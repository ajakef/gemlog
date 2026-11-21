"""
Checks Aspen logger config files (CONFIG.TXT, and config sequence files CONFIG01.TXT, CONFIG02.TXT, ...) before a
deployment, with the same rules the logger uses (Aspen firmware 1.0.0-34 and later). What each setting does:
config_file_aspen.md in the gemlog repository (https://github.com/ajakef/gemlog/blob/main/config_file_aspen.md).

From Python:
    import gemlog
    result = gemlog.check_aspen_config('E:/CONFIG.TXT')   # prints a report; result.ok is False if the logger would
                                                          # report a mistake (red LED: 6 flashes)
    gemlog.check_aspen_config('E:/')                       # a whole SD card: CONFIG.TXT and any sequence files
From a terminal:
    check_aspen_config E:/CONFIG.TXT

THIS MODULE MUST AGREE WITH THE FIRMWARE (aspen-firmware: src/firmware_logic.cpp, other_functions.h; see
"---- CONFIG.TXT ----" in src/firmware_logic.h). The settings table (names, ranges, defaults, words) is read from
gemlog/data/aspen_config_spec.json, a copy of the firmware's docs/config_spec.json, which the firmware's host tests
write from its own code together with example lines and files and what the firmware makes of them;
tests/test_aspen_config.py checks this module against those examples. After a firmware change to how config files are
read, run aspen-firmware/tools/sync_config_spec.sh and make the tests here pass.

Files are handled as bytes (each byte one character, as latin-1 reads it), so that line lengths and odd characters
count the way they do on the logger.
"""
import argparse
import json
import os
import re
import sys

_SPEC_PATH = os.path.join(os.path.dirname(__file__), 'data', 'aspen_config_spec.json')
_spec = None


def _load_spec():
    global _spec
    if _spec is None:
        with open(_SPEC_PATH, encoding = 'utf-8') as f:
            _spec = json.load(f)
        _spec['by_name'] = {s['name']: s for s in _spec['settings']}
    return _spec


_SPACES = ' \t\r\n'
_BOM = '\xef\xbb\xbf'
_NAME_MAX = 23   # ConfigFileLine.name holds 23 characters
_VALUE_MAX = 39  # ConfigFileLine.value holds 39
_FIRMWARE_RULES = '1.0.0-34'

# Positional (original) format: the order of the values, with names as the firmware's messages give them.
_POSITIONAL = ['channel 1 gain', 'channel 2 gain', 'channel 3 gain', 'channel 4 gain', 'sample rate', 'gps_mode',
               'gps_cycle_minutes', 'gps_quota', 'use_compression', 'network', 'station', 'location',
               'channel', 'channel', 'channel', 'channel', 'checkin_minutes', 'checkin_hour']
_POSITIONAL_LABELS = ['channel 1 gain', 'channel 2 gain', 'channel 3 gain', 'channel 4 gain', 'sample rate',
                      'gps_mode', 'gps_cycle_minutes', 'gps_quota', 'use_compression', 'network code',
                      'station code', 'location code', 'channel 1 code', 'channel 2 code', 'channel 3 code',
                      'channel 4 code', 'checkin_minutes', 'checkin_hour']
# Words in a positional file's comments that name each setting, to spot a missing or extra line.
_COMMENT_PATTERNS = ([re.compile(r'channel\s*%d\s*gain' % i) for i in range(1, 5)] +
                     [re.compile(p) for p in [r'sample[ _]rate', r'gps[ _]mode', r'cycle', r'quota', r'compression',
                                              r'network', r'station', r'location']] +
                     [re.compile(r'channel\s*%d\s*code' % i) for i in range(1, 5)] +
                     [re.compile(r'checkin_minutes|check-in.*minutes'), re.compile(r'checkin_hour|hour')])
_INT32_MIN = -2147483648


def _c_string(s):
    """The text as the firmware sees it in a C string: up to the first 0 byte."""
    i = s.find('\0')
    return s if i < 0 else s[:i]


def _ascii_lower(s):
    return ''.join(chr(ord(c) + 32) if 'A' <= c <= 'Z' else c for c in s)


def _is_digits(s):
    return s != '' and all('0' <= c <= '9' for c in s)


def _whole_number(text):
    """parseWholeNumber: an optional '-' and 1-9 digits, nothing else; None otherwise."""
    digits = text[1:] if text.startswith('-') else text
    if not _is_digits(digits) or len(digits) > 9:
        return None
    return -int(digits) if text.startswith('-') else int(digits)


def split_lines(text):
    """
    Splits a config file's text (one character per byte) into lines as the logger reads them (ConfigLineAssembler):
    a line ends at '\\n'; '\\r' is dropped anywhere; '#' starts a comment; only the first line_max characters before the
    comment are kept. Returns a list of (setting text, cut, line as written) with cut True if the line had more.
    """
    line_max = _load_spec()['line_max']
    lines = []
    kept, shown = [], []
    cut = in_comment = started = False
    for c in text:
        started = True
        if c == '\n':
            lines.append((_c_string(''.join(kept)), cut, ''.join(shown)))
            kept, shown = [], []
            cut = in_comment = started = False
            continue
        if c != '\r':
            shown.append(c)
        if c == '#':
            in_comment = True
        if c == '\r' or in_comment:
            continue
        if len(kept) < line_max:
            kept.append(c)
        else:
            cut = True
    if started:
        lines.append((_c_string(''.join(kept)), cut, ''.join(shown)))
    return lines


def detect_format(first_line):
    """detectConfigFormat: 'utf16', 'positional' (starts with a digit, '-' or '+'), or 'named'."""
    if first_line.startswith('\xff\xfe') or first_line.startswith('\xfe\xff'):
        return 'utf16'
    p = first_line[3:] if first_line.startswith(_BOM) else first_line
    p = p.lstrip(_SPACES)
    if p[:1].isascii() and (p[:1].isdigit() or p[:1] in ('-', '+')):
        return 'positional'
    return 'named'


def _name_distance(a, b):
    if len(a) > 24 or len(b) > 24:
        return 99
    row = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        diagonal, row[0] = row[0], i
        for j in range(1, len(b) + 1):
            above = row[j]
            row[j] = min(diagonal + (a[i - 1] != b[j - 1]), row[j - 1] + 1, above + 1)
            diagonal = above
    return row[len(b)]


def _closest_setting(name):
    test = name.startswith('test_')
    best, best_distance = None, 99
    for s in _load_spec()['settings']:
        candidate = s['name']
        if candidate.startswith('test_') != test:
            continue
        d = 1 if (len(name) >= 3 and candidate.startswith(name)) else _name_distance(name, candidate)
        if d < best_distance:
            best, best_distance = candidate, d
    return best if best_distance <= (3 if len(name) >= 8 else 2) else None


def _seed_code(text, max_len, exact):
    """parseSeedCode: letters and digits (made uppercase), up to max_len (exactly, if exact); None if not."""
    code = []
    for c in text:
        c = chr(ord(c) - 32) if 'a' <= c <= 'z' else c
        if not ('A' <= c <= 'Z' or '0' <= c <= '9') or len(code) >= max_len:
            return None
        code.append(c)
    if not code or (exact and len(code) != max_len):
        return None
    return ''.join(code)


def _number_value(setting, text):
    """A number setting's value (a whole number or one of its words), or None if it isn't allowed."""
    x = _whole_number(text)
    if x is None:
        x = setting.get('words', {}).get(_ascii_lower(text))
    if x is None or x < setting['min'] or x > setting['max']:
        return None
    if 'values' in setting and x not in setting['values']:
        return None
    return x


def _gain_list(text):
    """parseGainList: 1-4 comma-separated allowed gains (missing ones 0), not all 0; None if not."""
    allowed = _load_spec()['by_name']['gains']['values']
    fields = text.split(',')
    if len(fields) > 4:
        return None
    gains = []
    for f in fields:
        if not _is_digits(f):
            return None
        value = 0
        for c in f:  # the firmware gives up once a field passes 255
            value = value * 10 + int(c)
            if value > 255:
                return None
        if value not in allowed:
            return None
        gains.append(value)
    gains += [0] * (4 - len(gains))
    if not any(gains):
        return None
    return gains


def _code_assignment(name, text):
    """A SEED code setting (network, station, location, channel1-4): [name, code or '-'], or None."""
    s = _load_spec()['by_name'][name]
    if text == '-':
        return [name, '-']
    code = _seed_code(text, s['max_length'], s['exact'])
    return None if code is None else [name, code]


def parse_line(line):
    """
    What the logger makes of one line of a named config file (parseConfigFileLine), comment already removed. Returns
    a dict: status ('ok', 'blank', 'no_equals', 'unknown_name', 'test_setting', 'extra_text' or 'bad_value'), name
    and value (as written, for messages), assignments (a list of [setting, value]: gain1-4, channel1-4, network,
    station, location take '-' for the default, everything else a number), and suggestion (for unknown names).
    """
    spec = _load_spec()
    out = {'status': 'blank', 'name': '', 'value': '', 'assignments': [], 'suggestion': None}
    # a non-breaking space (bytes C2 A0) is a space; the firmware copies at most line_max characters
    chars, i = [], 0
    while i < len(line) and len(chars) < spec['line_max']:
        if line[i] == '\xc2' and line[i + 1:i + 2] == '\xa0':
            chars.append(' ')
            i += 2
        else:
            chars.append(line[i])
            i += 1
    text = ''.join(chars)
    p = text[3:] if text.startswith(_BOM) else text
    p = p.lstrip(_SPACES)
    if p == '':
        return out
    eq = p.find('=')
    out['name'] = _ascii_lower((p if eq < 0 else p[:eq]).rstrip(_SPACES)[:_NAME_MAX])
    name = out['name']
    if eq < 0:
        out['status'] = 'no_equals'
        return out
    out['value'] = p[eq + 1:].strip(_SPACES)[:_VALUE_MAX]
    setting = spec['by_name'].get(name)
    if setting is None:
        out['status'] = 'unknown_name'
        out['suggestion'] = _closest_setting(name)
        return out
    test = (not setting['config_txt']) and name.startswith('test_')

    # the value with spaces removed around commas; spaces anywhere else mean two values
    value, pending_space, extra = [], False, False
    for c in p[eq + 1:]:
        if c in _SPACES:
            pending_space = True
            continue
        if pending_space and value and value[-1] != ',' and c != ',':
            extra = True
        pending_space = False
        if len(value) >= 63:
            out['status'] = 'bad_value'
            return out
        value.append(c)
    value = ''.join(value)
    if value == '':
        out['status'] = 'bad_value'
        return out
    if extra:
        out['status'] = 'extra_text'
        return out

    kind = setting['type']
    items = None
    if value == '-':
        if kind == 'gains':
            items = [['gain%d' % (i + 1), g] for i, g in enumerate(setting['default'])]
        elif kind == 'channels':
            items = [['channel%d' % (i + 1), '-'] for i in range(4)]
        elif kind == 'code':
            items = [[name, '-']]
        elif name != 'run_minutes':  # run_minutes has no default
            items = [[name, setting['default']]]
    elif kind == 'gains':
        gains = _gain_list(value)
        if gains is not None:
            items = [['gain%d' % (i + 1), g] for i, g in enumerate(gains)]
    elif kind == 'channels':
        entries = value.split(',')
        if len(entries) <= 4 and all(len(e) <= 8 for e in entries):
            items = [_code_assignment('channel%d' % (i + 1), e) for i, e in enumerate(entries)]
            if None in items:
                items = None
    elif kind == 'code':
        item = _code_assignment(name, value)
        items = None if item is None else [item]
    elif name == 'run_minutes':
        x = _whole_number(value)
        if x is not None and setting['min'] <= x <= setting['max']:
            items = [[name, x]]
    else:
        x = _number_value(setting, value)
        if x is not None:
            items = [[name, x]]
    if items is None:
        out['status'] = 'bad_value'
        return out
    out['assignments'] = items
    out['status'] = 'test_setting' if test else 'ok'
    return out


def _setting_ids(assignments):
    return {a[0] for a in assignments}


def _problem_message(problem, parsed, file_name, line_number):
    """The logger's message for a problem (formatConfigProblem)."""
    spec = _load_spec()
    name, value = parsed['name'], parsed['value']
    allowed = spec['by_name'][name]['allowed'] if name in spec['by_name'] else ''
    if problem == 'no_equals':
        body = '"%s" has no \'=\' (expected name = value)' % name
    elif problem == 'unknown_name':
        body = '"%s" is not a setting' % name
        if parsed['suggestion']:
            body += ' (did you mean %s?)' % parsed['suggestion']
    elif problem == 'extra_text':
        body = '%s = %s has more than one value (allowed: %s)' % (name, value, allowed)
    elif problem == 'bad_value':
        body = "%s = %s isn't allowed (allowed: %s)" % (name, value, allowed)
    elif problem == 'test_setting_here':
        body = ('%s is a test setting: it works only in config sequence files (CONFIG01.TXT, ...) and CF: serial '
                'commands' % name)
    elif problem == 'run_minutes_here':
        body = 'run_minutes works only in config sequence files (CONFIG01.TXT, CONFIG02.TXT, ...)'
    elif problem == 'set_twice':
        body = '%s is set twice' % name
    else:  # too_long
        body = 'longer than %d characters before any comment' % spec['line_max']
    return 'Error: %s line %d: %s; %s' % (file_name, line_number, body,
                                         'this line is used' if problem == 'set_twice' else 'ignored')


def read_named(lines, sequence_file, file_name = 'CONFIG.TXT'):
    """
    Reads a named config file's lines (from split_lines) as NamedConfigReader does. Returns (problems, messages,
    settings, has_run_minutes): problems a list of [line number, kind]; messages the logger's "Error: ..." lines;
    settings the last value given for each setting (as in parse_line's assignments).
    """
    problems, messages, settings = [], [], {}
    seen = set()
    has_run_minutes = False
    for number, (text, cut, _) in enumerate(lines, start = 1):
        parsed = parse_line(text)
        problem = None
        apply = False
        info = _load_spec()['by_name'].get(parsed['name'])
        test = info is not None and not info['config_txt'] and parsed['name'].startswith('test_')
        status = parsed['status']
        if cut:
            problem = 'too_long'
        elif status == 'blank':
            pass
        elif status in ('no_equals', 'unknown_name'):
            problem = status
        elif status in ('extra_text', 'bad_value'):
            problem = 'test_setting_here' if (test and not sequence_file) else status
        elif status == 'test_setting' and not sequence_file:
            problem = 'test_setting_here'
        else:
            run_minutes = len(parsed['assignments']) == 1 and parsed['assignments'][0][0] == 'run_minutes'
            if run_minutes and not sequence_file:
                problem = 'run_minutes_here'
            else:
                has_run_minutes = has_run_minutes or run_minutes
                ids = _setting_ids(parsed['assignments'])
                if ids & seen:
                    problem = 'set_twice'
                seen |= ids
                apply = True
        if apply:
            for setting, value in parsed['assignments']:
                settings[setting] = value
        if problem is not None:
            problems.append([number, problem])
            messages.append(_problem_message(problem, parsed, file_name, number))
    return problems, messages, settings, has_run_minutes


def file_format(lines):
    """The file's format from its first line with anything on it (besides a comment), as the logger decides it."""
    for text, _, _ in lines:
        if parse_line(text)['status'] != 'blank':
            return detect_format(text)
    return 'named'


def sequence_number(name):
    """The number n in a config sequence file's name (CONFIG<n>.TXT, 1-3 digits, any case), or -1."""
    m = re.fullmatch(r'CONFIG([0-9]{1,3})\.TXT', name, flags = re.IGNORECASE | re.ASCII)
    return int(m.group(1)) if m else -1


def _strtol(text):
    """strtol as the positional reader uses it: the number at the start, or None if there is none ('missing')."""
    m = re.match(r'[+-]?[0-9]+', text)
    if m is None:
        return None
    x = max(min(int(m.group()), 2147483647), _INT32_MIN)
    return None if x == _INT32_MIN else x


def read_positional(text, file_name = 'CONFIG.TXT'):
    """
    Reads a positional (original format) file as ReadConfigPositional does. Returns (messages, warnings, settings,
    line_notes): the logger's error messages; warnings about lines the logger silently ignores; the settings in the
    same form as read_named's; and {line number: what the line is read as}.
    """
    spec = _load_spec()
    defaults = spec['by_name']
    values = []  # (line number, value text as the firmware sees it, comment)
    for number, raw in enumerate(text.split('\n'), start = 1):
        hash_at = raw.find('#')
        comment = raw[hash_at + 1:] if hash_at >= 0 else ''
        kept = ''.join(c for c in (raw if hash_at < 0 else raw[:hash_at]) if c not in ' \r\t')[:198]
        if kept:
            values.append((number, _c_string(kept), _ascii_lower(comment)))
    messages, warnings, settings, notes = [], [], {}, {}

    def invalid(x, used, name):
        messages.append('Error: invalid %s %d in %s; using %d' % (name, x, file_name, used))

    gains = list(defaults['gains']['default'])
    for index, (number, value, comment) in enumerate(values):
        if index >= len(_POSITIONAL):
            warnings.append('line %d: an extra line: the logger reads only the first %d values, so it ignores this'
                            % (number, len(_POSITIONAL)))
            continue
        label = _POSITIONAL_LABELS[index]
        notes[number] = label
        if comment.strip() and not _COMMENT_PATTERNS[index].search(comment):
            others = [k for k, pattern in enumerate(_COMMENT_PATTERNS) if pattern.search(comment)]
            if others:
                warnings.append('line %d: read as %s, but its comment says %s: is a line above it missing or extra?'
                                % (number, label, _POSITIONAL_LABELS[others[0]]))
        if 9 <= index <= 15:  # SEED codes
            if value == '-':
                continue
            name = ['network', 'station', 'location', 'channel1', 'channel2', 'channel3', 'channel4'][index - 9]
            s = defaults[name]
            code = _seed_code(value, s['max_length'], s['exact'])
            if code is None:
                messages.append('Error: invalid %s code "%s" in %s; using default' % (_POSITIONAL[index], value,
                                                                                      file_name))
            else:
                settings[name] = code
            continue
        x = _strtol(value)
        if x is None:
            warnings.append('line %d: "%s" isn\'t a number, so the logger ignores it and uses the default %s'
                            % (number, value, label))
            continue
        if not re.fullmatch(r'[+-]?[0-9]+', value):
            warnings.append('line %d: the logger reads only the number at the start of "%s"' % (number, value))
        if index < 4:
            if x in defaults['gains']['values']:
                gains[index] = x
            else:
                invalid(x, gains[index], _POSITIONAL[index])
        elif index == 4:
            if x in defaults['sample_rate']['values']:
                settings['sample_rate'] = x
            else:
                invalid(x, defaults['sample_rate']['default'], 'sample rate')
        elif index == 5:
            if 1 <= x <= 4:
                settings['gps_mode'] = x
            else:
                invalid(x, defaults['gps_mode']['default'], 'gps_mode')
        elif index in (6, 7):
            name = 'gps_cycle_minutes' if index == 6 else 'gps_quota'
            if 0 < x <= 255:
                settings[name] = x
            else:
                invalid(x, defaults[name]['default'], name)
        elif index == 8:
            if x not in (0, 1):
                invalid(x, 1, 'use_compression')  # normalizeUseCompression: invalid means 1 (it changes nothing now)
        elif index == 16:
            if 0 <= x <= 30:
                settings['checkin_minutes'] = x
            else:
                invalid(x, defaults['checkin_minutes']['default'], 'checkin_minutes')
        elif index == 17:
            if 0 <= x <= 23:
                settings['checkin_hour'] = x
            else:
                invalid(x, defaults['checkin_hour']['default'], 'checkin_hour')
    if not any(gains):
        messages.append('Error: all channels set to gain 0 in %s; using default gains' % file_name)
        gains = list(defaults['gains']['default'])
    for i, g in enumerate(gains):
        settings['gain%d' % (i + 1)] = g
    return messages, warnings, settings, notes


def _default_channel_code(channel, sample_rate):
    """defaultChannelCode: HDF, EHZ, EHN, EHE at 100-200 sps; band C/D at 400-800; F/G at 1000 and up."""
    if sample_rate >= 1000:
        infrasound, seismic = 'F', 'G'
    elif sample_rate >= 250:
        infrasound, seismic = 'C', 'D'
    else:
        infrasound, seismic = 'H', 'E'
    return infrasound + 'DF' if channel == 0 else seismic + 'H' + 'ZNE'[channel - 1]


class AspenConfigCheck:
    """
    The result of check_aspen_config for one file.

    Attributes
    ----------
    file_name : the file's name (CONFIG.TXT, CONFIG01.TXT, ...)
    format : 'named', 'positional', or 'utf16' (unreadable: every setting takes its default)
    errors : the logger's own messages for what it would report ("Error: ..."); any of them makes the red LED
        flash 6 times, and the line (or value) concerned is ignored
    warnings : things the logger accepts without a message that are probably not what was meant
    notes : other things worth knowing
    settings : every setting as the logger would use it (defaults filled in)
    lines : (line number, line as written, what it is: a setting, a problem, or '') for showing the file
    sequence : for a config sequence file, whether the logger would run it (True/False); None otherwise
    ok : True if there are no errors
    """

    def __init__(self, file_name):
        self.file_name = file_name
        self.format = 'named'
        self.errors = []
        self.warnings = []
        self.notes = []
        self.settings = {}
        self.lines = []
        self.sequence = None

    @property
    def ok(self):
        return not self.errors

    def report(self):
        out = ['%s: %s format' % (self.file_name, {'named': 'named (name = value)', 'positional':
                                                  'positional (original)', 'utf16': 'UTF-16, unreadable'}[self.format])]
        if self.sequence is not None:
            out.append('  Config sequence file: ' + ('the logger runs it' if self.sequence else
                                                     'the logger IGNORES it (no valid run_minutes line)'))
        for title, items in (('Errors (the logger reports these: red LED 6 flashes)', self.errors),
                             ('Warnings', self.warnings), ('Notes', self.notes)):
            if items:
                out.append('  %s:' % title)
                out += ['    ' + item for item in items]
        if self.settings and self.sequence is not False:
            out.append('  Settings the logger would use:')
            s = self.settings
            out.append('    gains %s; %s sps; gps_mode %s, cycle %s min, quota %s, settle %s s' % (
                ', '.join(str(s['gain%d' % i]) for i in range(1, 5)), s['sample_rate'], s['gps_mode'],
                s['gps_cycle_minutes'], s['gps_quota'], s['gps_settle_s']))
            out.append('    SEED codes %s.%s.%s, channels %s' % (
                s['network'] or '(blank)', s['station'], s['location'] or '(blank)',
                ', '.join(s['channel%d' % i] for i in range(1, 5))))
            out.append('    files %s min; H lines every %s s; serial data/gps/health %s/%s/%s; bluetooth %s; '
                       'check-in %s min at %s h UTC' % (
                           s['file_minutes'], s['health_s'], s['serial_data'], s['serial_gps'], s['serial_health'],
                           s['bluetooth'], s['checkin_minutes'], s['checkin_hour']))
            extra = ['%s %s' % (k, v) for k, v in sorted(s.items()) if k.startswith('test_') or k == 'run_minutes']
            if extra:
                out.append('    ' + '; '.join(extra))
        out.append('  ' + ('OK: no errors' if self.ok else '%d error(s)' % len(self.errors)))
        return '\n'.join(out)

    def __str__(self):
        return self.report()

    def __repr__(self):
        return '<AspenConfigCheck %s: %s, %d errors, %d warnings>' % (self.file_name, self.format, len(self.errors),
                                                                      len(self.warnings))


def _check_bytes(data, file_name):
    spec = _load_spec()
    result = AspenConfigCheck(file_name)
    text = data.decode('latin-1')
    lines = split_lines(text)
    result.format = file_format(lines)
    is_sequence = sequence_number(file_name) >= 0
    settings = {}
    line_notes = {}
    if result.format == 'utf16':
        result.errors.append('Error: %s is saved as UTF-16 ("Unicode"), which the logger can\'t read: save it as '
                             'plain text (UTF-8 or ANSI); using the defaults' % file_name)
    elif result.format == 'positional':
        messages, warnings, settings, line_notes = read_positional(text, file_name)
        result.errors += messages
        result.warnings += warnings
        result.notes.append('This is the original positional format (one value per line, in a fixed order). The '
                            'named format (name = value) is easier to edit: see config_file_aspen.md.')
    else:
        name_in_messages = file_name if is_sequence else 'CONFIG.TXT'
        problems, messages, settings, has_run_minutes = read_named(lines, is_sequence, name_in_messages)
        result.errors += messages
        for number, problem in problems:
            line_notes[number] = 'problem'
        if is_sequence:
            result.sequence = has_run_minutes
        result.notes.append('Named files need Aspen firmware 1.0.0-30 or later (older firmware silently uses the '
                            'defaults); these checks follow firmware %s and later.' % _FIRMWARE_RULES)
    if is_sequence and result.sequence is None:
        result.sequence = False
    if is_sequence and not result.sequence:
        result.warnings.append('The logger ignores this file: a config sequence file needs a valid run_minutes line '
                               'in the named format (Aspen firmware %s and later).' % _FIRMWARE_RULES)
    elif is_sequence:
        result.notes.append('Config sequence file: the logger uses it instead of CONFIG.TXT for %s min, then renames '
                            'it %s and restarts with the next one.' % (settings.get('run_minutes'),
                                                                       re.sub(r'(?i)\.txt$', '_DONE.TXT', file_name)))
    for number, (setting_text, cut, shown) in enumerate(lines, start = 1):
        result.lines.append((number, shown, line_notes.get(number, '')))

    # The settings as the logger would use them: defaults, then the file, then what setup() does at startup.
    effective = {}
    for s in spec['settings']:
        if s['type'] == 'number' and 'default' in s and s['config_txt']:
            effective[s['name']] = s['default']
    for i, g in enumerate(spec['by_name']['gains']['default']):
        effective['gain%d' % (i + 1)] = g
    for name in ('network', 'station', 'location', 'channel1', 'channel2', 'channel3', 'channel4'):
        effective[name] = '-'
    effective.update(settings)
    gains = [effective['gain%d' % i] for i in range(1, 5)]
    rate = effective['sample_rate']
    if rate >= 4000 and (gains[0] == 0 or any(gains[1:])):
        result.errors.append('Error: %d sps records channel 1 only; channels 2-4 off' % rate)
        gains = [gains[0] or 128, 0, 0, 0]
        for i, g in enumerate(gains):
            effective['gain%d' % (i + 1)] = g
    for name in ('network', 'location'):
        if effective[name] == '-':
            effective[name] = ''
    if effective['station'] == '-':
        effective['station'] = "(the logger's serial number)"
    for i in range(4):
        if effective['channel%d' % (i + 1)] == '-':
            effective['channel%d' % (i + 1)] = _default_channel_code(i, rate)
    codes = [effective['channel%d' % (i + 1)] for i in range(4)]
    for i in range(1, 4):
        for j in range(i):
            if gains[i] and gains[j] and codes[i] == codes[j]:
                result.errors.append('Error: channel %d has the same code as channel %d; using XX%d' % (i + 1, j + 1,
                                                                                                       i + 1))
                codes[i] = 'XX%d' % (i + 1)
                effective['channel%d' % (i + 1)] = codes[i]
                break
    if rate >= 4000 and effective['bluetooth']:
        result.notes.append('Bluetooth stays off at 4000 sps.')
    result.settings = effective
    return result


def check_aspen_config(path = None, text = None, name = None, verbose = True):
    """
    Check an Aspen config file (or every config file on an SD card) with the logger's own rules.

    Parameters
    ----------
    path : str
        A config file (CONFIG.TXT, CONFIG01.TXT, ...), or a folder such as an SD card's root, whose CONFIG.TXT and
        config sequence files are all checked.
    text : str or bytes
        Instead of path: the file's contents (a str is taken as typed, i.e. saved as UTF-8).
    name : str
        The file's name when text is given (default 'CONFIG.TXT'); a name like CONFIG01.TXT makes it a config
        sequence file.
    verbose : bool
        Print the report (default True).

    Returns
    -------
    AspenConfigCheck for one file (see its attributes: ok, errors, warnings, notes, settings, ...), or a list of them
    for a folder (in the order the logger would use them).
    """
    if text is not None:
        data = text.encode('utf-8') if isinstance(text, str) else bytes(text)
        result = _check_bytes(data, name or 'CONFIG.TXT')
        results = [result]
    elif path is not None and os.path.isdir(path):
        files = sorted(os.listdir(path))
        sequence = sorted((sequence_number(f), f) for f in files if sequence_number(f) >= 0)
        main = [f for f in files if f.upper() == 'CONFIG.TXT']
        results = []
        for _, f in sequence + [(None, m) for m in main]:
            with open(os.path.join(path, f), 'rb') as handle:
                results.append(_check_bytes(handle.read(), f))
        used = next((r for r in results if r.sequence), None)
        if used is None:
            used = next((r for r in results if r.file_name.upper() == 'CONFIG.TXT'), None)
        summary = ('The logger would use %s first.' % used.file_name if used is not None else
                   'No CONFIG.TXT and no config sequence file: the logger would use its defaults.')
        if 'CONFIG.NEW' in [f.upper() for f in files]:
            summary += (' CONFIG.NEW is left from a save (SVC) that was cut off; the logger finishes it at startup if '
                        'there is no CONFIG.TXT.')
        if verbose:
            print('\n\n'.join(r.report() for r in results) + ('\n\n' if results else '') + summary)
        return results
    elif path is not None:
        with open(path, 'rb') as handle:
            result = _check_bytes(handle.read(), name or os.path.basename(path))
        results = [result]
    else:
        raise ValueError('give path or text')
    if verbose:
        print(results[0].report())
    return results[0]


def main(argv = None):
    parser = argparse.ArgumentParser(
        description = 'Check Aspen logger config files (CONFIG.TXT, CONFIG01.TXT, ...) with the logger\'s own rules. '
                      'Exit status 1 if the logger would report a mistake in any of them.')
    parser.add_argument('paths', nargs = '+', help = 'config files, or folders such as an SD card\'s root')
    args = parser.parse_args(argv)
    ok = True
    for i, path in enumerate(args.paths):
        if i:
            print()
        result = check_aspen_config(path)
        for r in (result if isinstance(result, list) else [result]):
            ok = ok and r.ok
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())

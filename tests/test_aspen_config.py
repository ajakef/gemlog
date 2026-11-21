"""
Tests for gemlog.aspen_config, the Aspen config file checker. It must agree with the logger's firmware: the expected
results below come from gemlog/data/aspen_config_spec.json, which the firmware's host tests write from its own code
(aspen-firmware docs/config_spec.json, copied here by aspen-firmware/tools/sync_config_spec.sh).
"""
import json
import os
import pytest
from gemlog import aspen_config
from gemlog.aspen_config import check_aspen_config, parse_line, split_lines, file_format, read_named, sequence_number

SPEC_PATH = os.path.join(os.path.dirname(aspen_config.__file__), 'data', 'aspen_config_spec.json')
with open(SPEC_PATH, encoding = 'utf-8') as f:
    SPEC = json.load(f)
FIRMWARE_SPEC = os.path.join(os.path.dirname(__file__), '..', '..', 'aspen-firmware', 'docs', 'config_spec.json')


@pytest.mark.skipif(not os.path.exists(FIRMWARE_SPEC), reason = 'aspen-firmware is not next to gemlog')
def test_spec_copy_matches_the_firmware():
    with open(FIRMWARE_SPEC, encoding = 'utf-8') as f:
        assert json.load(f) == SPEC, 'run aspen-firmware/tools/sync_config_spec.sh'


@pytest.mark.parametrize('case', SPEC['line_cases'], ids = lambda c: repr(c['line'])[:40])
def test_lines_read_as_the_firmware_reads_them(case):
    parsed = parse_line(case['line'])
    assert parsed['status'] == case['status']
    if 'assignments' in case:
        assert parsed['assignments'] == case['assignments']
    if 'suggestion' in case:
        assert parsed['suggestion'] == case['suggestion']


@pytest.mark.parametrize('case', SPEC['file_cases'], ids = lambda c: c['name'] + ' ' + repr(c['text'])[:30])
def test_files_read_as_the_firmware_reads_them(case):
    lines = split_lines(case['text'])
    assert file_format(lines) == case['format']
    if case['format'] != 'named':
        return
    sequence = sequence_number(case['name']) >= 0
    problems, messages, settings, has_run_minutes = read_named(lines, sequence,
                                                               'CONFIG01.TXT' if sequence else 'CONFIG.TXT')
    assert problems == case['problems']
    assert messages == case['messages']
    assert settings == case['settings']
    if 'has_run_minutes' in case:
        assert has_run_minutes == case['has_run_minutes']


def test_standard_config_file_is_the_firmwares_and_clean():
    path = os.path.join(os.path.dirname(aspen_config.__file__), 'data', 'aspen_CONFIG.TXT')
    with open(path, 'rb') as f:
        data = f.read()
    assert data.decode('latin-1') == SPEC['file_cases'][0]['text']  # the firmware repository's CONFIG.TXT
    assert check_aspen_config(path, verbose = False).ok


def test_sequence_file_names():
    for name, number in SPEC['sequence_names']:
        assert sequence_number(name) == number, name


def test_repository_config_file_is_clean():
    result = check_aspen_config(text = SPEC['file_cases'][0]['text'].encode('latin-1'), verbose = False)
    assert result.format == 'named'
    assert result.ok
    assert result.settings['sample_rate'] == 200
    assert [result.settings['gain%d' % i] for i in range(1, 5)] == [128, 128, 128, 128]
    assert result.settings['channel1'] == 'HDF'
    assert 'serial number' in result.settings['station']


def test_named_file_mistakes_are_the_loggers_messages():
    result = check_aspen_config(text = 'sample_rate = 300\ngps_mod = 1\nbluetooth = off\n', verbose = False)
    assert result.errors == [
        "Error: CONFIG.TXT line 1: sample_rate = 300 isn't allowed (allowed: 100, 200, 400, 800, 1000, 2000 or 4000; "
        "default 200); ignored",
        'Error: CONFIG.TXT line 2: "gps_mod" is not a setting (did you mean gps_mode?); ignored']
    assert result.settings['sample_rate'] == 200 and result.settings['bluetooth'] == 0
    assert not result.ok


def test_startup_checks():
    result = check_aspen_config(text = 'sample_rate = 4000\ngains = 128, 32\n', verbose = False)
    assert 'Error: 4000 sps records channel 1 only; channels 2-4 off' in result.errors
    assert result.settings['gain2'] == 0
    result = check_aspen_config(text = 'channels = HDF, HDF\n', verbose = False)
    assert 'Error: channel 2 has the same code as channel 1; using XX2' in result.errors
    result = check_aspen_config(text = 'sample_rate = 400\n', verbose = False)
    assert [result.settings['channel%d' % i] for i in range(1, 5)] == ['CDF', 'DHZ', 'DHN', 'DHE']


def test_positional_files():
    good = SPEC['file_cases'][1]['text'].encode('latin-1')  # the firmware repository's positional example
    result = check_aspen_config(text = good, verbose = False)
    assert result.format == 'positional'
    assert result.ok, result.errors
    assert result.settings['sample_rate'] == 400 and result.settings['gain2'] == 64
    bad = b'128 # channel 1 gain\n100 # channel 2 gain\n64\n64\n450\n1\n30\n20\n1\nXYZ\n'
    result = check_aspen_config(text = bad, verbose = False)
    assert result.errors == ['Error: invalid channel 2 gain 100 in CONFIG.TXT; using 128',
                             'Error: invalid sample rate 450 in CONFIG.TXT; using 200',
                             'Error: invalid network code "XYZ" in CONFIG.TXT; using default']
    # a deleted line shifts the rest: the comments give it away; a non-number is ignored without a message
    shifted = b'128 # channel 1 gain\n64 # channel 2 gain\n64 # channel 3 gain\n400 # sample rate\n1 # gps_mode\n'
    result = check_aspen_config(text = shifted, verbose = False)
    assert any('comment says sample rate' in w for w in result.warnings)
    result = check_aspen_config(text = b'128\nsample_rate = 400\n', verbose = False)
    assert result.ok and any("isn't a number" in w for w in result.warnings)


def test_utf16_and_first_line_mistakes():
    result = check_aspen_config(text = 'gains = 128\n'.encode('utf-16'), verbose = False)
    assert result.format == 'utf16' and not result.ok
    result = check_aspen_config(text = 'gains: 128\nsample_rate = 400\n', verbose = False)
    assert result.format == 'named'  # not read as a positional file
    assert result.settings['sample_rate'] == 400
    assert result.errors == ['Error: CONFIG.TXT line 1: "gains: 128" has no \'=\' (expected name = value); ignored']


def test_sequence_files(tmp_path):
    (tmp_path / 'CONFIG.TXT').write_text('sample_rate = 100\n')
    (tmp_path / 'CONFIG1.TXT').write_text('sample_rate = 400\n')  # e.g. an old copy: no run_minutes
    (tmp_path / 'CONFIG02.TXT').write_text('run_minutes = 60\ntest_card_full_s = 600\n')
    results = check_aspen_config(str(tmp_path), verbose = False)
    assert [r.file_name for r in results] == ['CONFIG1.TXT', 'CONFIG02.TXT', 'CONFIG.TXT']
    assert results[0].sequence is False and any('ignores this file' in w for w in results[0].warnings)
    assert results[1].sequence is True and results[1].ok
    assert results[1].settings['test_card_full_s'] == 600
    assert results[2].ok


def test_command_line(tmp_path, capsys):
    good = tmp_path / 'CONFIG.TXT'
    good.write_text('gps_mode = adaptive\n')
    assert aspen_config.main([str(good)]) == 0
    bad = tmp_path / 'bad' / 'CONFIG.TXT'
    bad.parent.mkdir()
    bad.write_text('gps_mode = 7\n')
    assert aspen_config.main([str(good), str(bad)]) == 1
    assert "gps_mode = 7 isn't allowed" in capsys.readouterr().out

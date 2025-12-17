"""
Custom exceptions for gemlog package.
"""


class EmptyRawFile(Exception):
    """Raised when the input file is empty"""
    pass


class CorruptRawFile(Exception):
    """Raised when the input file cannot be read for some reason"""
    pass


class CorruptRawFileNoGPS(CorruptRawFile):
    """Raised when the input file does not contain any GPS data"""
    pass


class CorruptRawFileInadequateGPS(CorruptRawFile):
    """Raised when the input file does not contain adequate GPS data"""
    pass


class CorruptRawFileDiscontinuousGPS(CorruptRawFile):
    """GPS timing discontinuity found in raw file due to receiving bad data from GPS; raw file must be repaired manually"""
    pass


class MissingRawFiles(Exception):
    """Raised when no input files are readable"""
    pass

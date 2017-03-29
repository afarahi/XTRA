import warnings

try:
    ImportWarning
except NameError:
    class ImportWarning(Warning):
        pass

try:
    from XTRA_Xray_SB_Class import Surface_Brightness_Class
except ImportError:
    warnings.warn("Warning: Cannot import XTRA_Xray_SB_Class",
                  category=ImportWarning)




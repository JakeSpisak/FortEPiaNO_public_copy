# wrapper for the main FortEPiaNO code
import six
from fortepianoWrapper import fortepianowrapper as fpw


def getVersion():
    """Print the version of the fortran code"""
    return fpw.getversion() if six.PY2 else str(fpw.getversion(), "utf-8")


if __name__ == "__main__":
    print("Wrapper for FortEPiaNO v%s" % getVersion())

from fortepianoWrapper import fortepianowrapper as fp

try:
    ver = str(fp.get_version(), "utf-8")
except TypeError:
    ver = fp.get_version()
print("--->Wrapper for FortEPiaNO v%s compiled and imported correctly!" % ver)

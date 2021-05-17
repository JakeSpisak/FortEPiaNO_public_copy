from fortepianoWrapper import fortepianowrapper as fp

try:
    ver = str(fp.getversion(), "utf-8")
except TypeError:
    ver = fp.getversion()
print("--->Wrapper for FortEPiaNO v%s compiled and imported correctly!" % ver)

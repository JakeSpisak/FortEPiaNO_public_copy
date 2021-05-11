import six
from fortepianoWrapper import fortepianowrapper as fp

print(
    "--->Wrapper for FortEPiaNO v%s compiled and imported correctly!"
    % (fp.getversion() if six.PY2 else str(fp.getversion(), "utf-8"))
)

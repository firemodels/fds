
import json
import urllib.request

FDSURL = "https://api.github.com/repos/firemodels/fds/releases"
SMVURL = "https://api.github.com/repos/firemodels/smv/releases"
BUNDLE = "FDS-6.10.1_SMV-6.10.1"
SMV = "SMV-6.10.5"
LINUX = "_lnx.sh"
OSX = "_osx.sh"
WIN = "_win.exe"


def get_downloads(filename: str, releases_url: str) -> int:
    """
    Retrieve the download count for a given asset filename from a GitHub releases feed.
    Returns 0 if the asset cannot be found.
    """
    request = urllib.request.Request(
        releases_url,
        headers={
            "Accept": "application/vnd.github+json",
            "User-Agent": "python-urllib"
        },
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        releases = json.load(response)

    for release in releases:
        for asset in release.get("assets", []):
            if asset.get("name") == filename:
                return int(asset.get("download_count", 0))
    return 0


print("bundle downloads")
print(f"FILE={BUNDLE}{LINUX} downloads={get_downloads(BUNDLE + LINUX, FDSURL)}")
print(f"FILE={BUNDLE}{OSX}   downloads={get_downloads(BUNDLE + OSX, FDSURL)}")
print(f"FILE={BUNDLE}{WIN}   downloads={get_downloads(BUNDLE + WIN, FDSURL)}")

print("")
print("smokeview downloads")
print(f"FILE={SMV}{LINUX} downloads={get_downloads(SMV + LINUX, SMVURL)}")
print(f"FILE={SMV}{OSX}   downloads={get_downloads(SMV + OSX, SMVURL)}")
print(f"FILE={SMV}{WIN}   downloads={get_downloads(SMV + WIN, SMVURL)}")


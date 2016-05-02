
if [ $1 == "mingw" ]
then
    EXT=.exe
    BUILD_TARGET=mingw_win_64
    PLAT=mingw
elif [ $1 == "gnu_linux" ]
then
    EXT=
    BUILD_TARGET=gnu_linux_64
    PLAT=linux
else
	exit 1
fi

SMV=../../SMV/Build/$BUILD_TARGET/smokeview_${PLAT}_64${EXT}
SMVOPTS=-killscript

$SMV $SMVOPTS -luascript test1.lua room_fire

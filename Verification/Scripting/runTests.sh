
if [ $1 == "mingw" ]
then
    EXT=.exe
    BUILD_TARGET=mingw_win_64
    PLAT=mingw
elif [ $1 == "gnu_linx" ]
then
    EXT=
    BUILD_TARGET=gnu_linx_64
    PLAT=linux
fi

SMV=../../SMV/Build/$BUILD_TARGET/smokeview_${PLAT}_64${EXT}
SMVOPTS=-killscript

$SMV $SMVOPTS -luascript test1.lua room_fire

#!/usr/bin/env bash

black='\E[30m'
red='\E[31m'
green='\E[32m'
yellow='\E[33m'
blue='\E[34m'
magenta='\E[35m'
cyan='\E[36m'
white='\E[37m'


cecho ()                     # Color-echo.
                             # Argument $1 = message
                             # Argument $2 = color
{
local default_msg="No message passed."
                             # Doesn't really need to be a local variable.

message=${1:-$default_msg}   # Defaults to default message.
color=${2:-$black}           # Defaults to black, if not specified.

  echo -en "$color"
  echo "$message"
  echo -en '\E[0m'                    # Reset to normal.

  return
}

# setup a basic test format
function runTest () {
    local expected=$1
    local test_name=$2
    local test_command=$3
    local test_outlog="$test_name.log"
    local test_outerr="$test_name.err"
    echo -n $test_name:
    $test_command > "$test_outlog" 2> "$test_outerr"
    exit_status=$?
    if [ $exit_status -eq $expected ]; then
        cecho "[OK]" $green
        rm "$test_outlog" "$test_outerr"
    else
        cecho "[Failed]" $red
        cat "$test_outerr"
        if [ $expected -eq 1 ]; then
            echo -n "   Smokeview should have returned an error with this input"
            echo " but did not"
        fi
    fi
    return $exit_status
}

# determine platform
UNAME=$(uname)
WORD_SIZE=64
# TODO: adjust for compiler
if [ "$UNAME" = "MINGW64_NT-10.0" ]
then
  PLATFORM=mingw
  BUILD_TARGET=mingw_win
  EXT=.exe
elif [ "$UNAME" = "Linux" ]
then
  PLATFORM=linux
  BUILD_TARGET=gnu_linux
else
    echo "platform could not be determined"
    exit 1
fi
echo PLATFORM=$PLATFORM
echo BUILD_TARGET=$BUILD_TARGET

if [ "$USE_SYSTEM_SMV" = "true" ]
then
    SMV=smokeview
    SMV=$(which "$SMV")
    SMV=$(readlink -f "$SMV") # use absolute path
else
    SMV=../../SMV/Build/smokeview/${BUILD_TARGET}_${WORD_SIZE}/smokeview_${PLATFORM}_${WORD_SIZE}${EXT}
    SMV=$(which "$SMV")
    SMV=$(readlink -f "$SMV") # use absolute path
fi

if [ "$USE_SYSTEM_FDS" = "true" ]
then
    FDS=fds
    FDS=$(which "$FDS")
    FDS=$(readlink -f "$FDS") # use absolute path
else
    FDS=fds
    # TODO: complete this
    #FDS=../../FDS_Compilation/smokeview/${BUILD_TARGET}_${WORD_SIZE}/smokeview_${PLATFORM}_${WORD_SIZE}${EXT}
    FDS=$(which "$FDS")
    FDS=$(readlink -f "$FDS") # use absolute path
fi

# simply run smokeview to determine if it is present
# if it is not present, abort testing
echo -n "running smokeview... "
"$SMV" -version > /dev/null
if [ $? -ne 0 ]
then
    echo "Smokeview could not be found at: $SMV"
    exit 1
else
    echo "present"
fi

# simply run fds to determine if it is present
# if it is not present, abort testing
echo -n "running fds... "
"$FDS" << EOF 2> /dev/null
EOF
if [ $? -ne 0 ]
then
    echo "FDS could not be found at: $FDS"
    exit 1
else
    echo present
fi

echo "building test data"
if [ -f test_outputs/room_fire/room_fire.stop ]
then
    echo "input data already exists"
else
    mkdir -p test_outputs/room_fire
    cp test_inputs/room_fire.fds test_outputs/room_fire
    (cd test_outputs/room_fire && "$FDS" room_fire.fds && touch room_fire.stop)
fi
echo "running test1.log"
(cd test_outputs/room_fire \
    && "$SMV" -killscript -luascript ../../tests/test1.lua \
        room_fire > test1.log)

echo "testing paths"
(mkdir -p test_outputs/path_testing/obs1/obs2/obs3/obs4/room_fire \
    && cp test_outputs/room_fire/room_fire* test_outputs/path_testing/obs1/obs2/obs3/obs4/room_fire \
    && cd test_outputs/path_testing/obs1/obs2/obs3 \
    && pwd \
    && "$SMV" -killscript -luascript ../../../../../tests/test1.lua \
        obs4/room_fire/room_fire.smv > test1.log)
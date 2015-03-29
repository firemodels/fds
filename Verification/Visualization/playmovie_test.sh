#!/bin/bash
movietype=.avi
if [ "`uname`" == "Darwin" ]; then
  movietype=.mp4
fi

ffplay movie_test$movietype
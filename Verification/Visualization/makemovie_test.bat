@echo off
erase movie_test*.png 2>Nul
smokeview -runscript movie_test 
ffmpeg -y -i movie_test_%%04d.png movie_test.avi
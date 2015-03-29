@echo off
erase movie_test*.png 2>Nul
echo making images
smokeview -runscript movie_test > Nul 2>Nul
echo making movie
ffmpeg -y -i movie_test_%%04d.png movie_test.avi
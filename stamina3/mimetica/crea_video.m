system("ffmpeg -framerate 1 -pattern_type glob -i './temp_img/Cell/*.png' -vf scale=1080:-1 video.mp4")
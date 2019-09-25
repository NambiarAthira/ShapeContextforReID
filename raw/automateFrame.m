
% Change the current folder to the folder of this m-file.
% (The line of code below is from Brett Shoelson of The Mathworks.)
% cd(fileparts(which(mfilename)));

% Define some top-level folder.
start_path = cd
topLevelFolder = uigetdir(start_path)

cd(topLevelFolder)
cd 1_1
copyfile('..\1_1\rgbframe*.png','C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\ATHIRA\HSVDatabase\1')

cd(topLevelFolder)
cd 1_2
copyfile('..\1_2\rgbframe*.png','C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\ATHIRA\HSVDatabase\1')

cd(topLevelFolder)
cd 1_3
copyfile('..\1_3\rgbframe*.png','C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\ATHIRA\HSVDatabase\1')

cd(topLevelFolder)
cd 1_4
copyfile('..\1_4\rgbframe*.png','C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\ATHIRA\HSVDatabase\1')

cd(topLevelFolder)
cd 1_5
copyfile('..\1_5\rgbframe*.png','C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\ATHIRA\HSVDatabase\1')

cd('C:\Users\AthiraNambiar\Desktop\PHD\codes\MoG\ATHIRA\HSVDatabase\1')
length(dir('rgbframe*.png'))

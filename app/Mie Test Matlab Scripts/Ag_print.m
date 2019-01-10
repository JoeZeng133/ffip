function Ag_print( fileID )
%AG_PRINT Summary of this function goes here
%   Detailed explanation goes here
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 3\n", 1.4447, 0, 1, 0);
fprintf(fileID, "{ Drude %e %e }\n", 1.3280e16/(2*pi), 9.1269e13/(2*pi));
fprintf(fileID, "{ CP %e %e %e %e}\n", -1.5951, 3.1288, 8.2749e15/(2*pi), 5.177e15/(2*pi));
fprintf(fileID, "{ CP %e %e %e %e}\n", 0.25261, -1.5066,  6.1998e15/(2*pi), 5.4126e14/(2*pi));
fprintf(fileID, "}\n");

end


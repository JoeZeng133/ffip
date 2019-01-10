function Au_print(fileID )
%AU_PRINT Summary of this function goes here
%   Detailed explanation goes here
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 3\n", 1.1431, 0, 1, 0);
fprintf(fileID, "{ Drude %e %e }\n", 1.3202e16/(2*pi), 1.0805e14/(2*pi));
fprintf(fileID, "{ CP %e %e %e %e}\n", 0.26698, -1.2371, 3.8711e15/(2*pi), 4.4642e14/(2*pi));
fprintf(fileID, "{ CP %e %e %e %e}\n", 3.0834, -1.0968, 4.1684e15/(2*pi), 2.3555e15/(2*pi));
fprintf(fileID, "}\n");
end


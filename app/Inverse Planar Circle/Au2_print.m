function Au2_print(fileID )
%AU_PRINT Summary of this function goes here
%   Detailed explanation goes here
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 1\n", 7.1431, 0, 1, 0);
fprintf(fileID, "{ Drude %e %e }\n", 1.3202e16/(2*pi), 1.0805e14/(2*pi));
fprintf(fileID, "}\n");
end


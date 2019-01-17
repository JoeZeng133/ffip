function Au2_print(fileID )
%AU_PRINT Summary of this function goes here
%   Detailed explanation goes here
fprintf(fileID, "{ ");
fprintf(fileID, "%e %e %e %e 2\n", 5.9673, 0, 1, 0);
fprintf(fileID, "{ Drude %e %e }\n", 2113.6e12, 15.92e12);
fprintf(fileID, '{ Lorentz %e %e %e}\n', 1.09, 650.07e12, 104.86e12);
fprintf(fileID, "}\n");
end


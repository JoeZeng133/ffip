function res = call_exe(func_name)
% call exe from utility or the workign directory 
% append utility path to PATH if not already 
full = mfilename('fullpath');
folder = fileparts(full);
path = getenv('PATH');
if ~contains(path, folder)
    setenv('PATH', [getenv('PATH') ';' folder]);
end

str = computer;
if contains(str, 'WIN')
   res = system([func_name, '.exe']);
else
   res = system(['./', func_name]);
end
end
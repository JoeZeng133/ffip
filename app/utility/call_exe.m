function res = call_exe(func_name)
    str = computer;
    if contains(str, 'WIN')
       res = system([func_name, '.exe']);
    else
       res = system(['./', func_name]);
    end
end
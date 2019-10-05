# creates formatoptions.rst
from ase.io.formats import all_formats, format2modulename
from ase.utils import import_module
with open('formatoptions.rst', 'w') as fd:
    print('============================================', file=fd)
    print('Format Specific Options', file=fd)
    print('============================================', file=fd)
    modules = [ format2modulename.get(form, form.replace('-','_')) for form in all_formats ]
    modules = list(set(modules))
    for module_name in sorted(modules):
        module = import_module('ase.io.'+module_name)
        read = None
        write = None
        if hasattr(module, 'read_'+module_name.replace('-','_')):
            read = 'read_'+module_name.replace('-','_')
        if hasattr(module, 'write_'+module_name.replace('-','_')):
            write = 'write_'+module_name.replace('-','_')
        print(module_name, file=fd)
        print('----------------------------------------', file=fd)
        if read:
            print('\n .. autofunction:: ase.io.{:}.{:}\n\n'.format(module_name,read), file=fd)
        if write:
            print('\n .. autofunction:: ase.io.{:}.{:}\n\n'.format(module_name,write), file=fd)
        if (not read) and (not write):
            print('\n No automatic documentation of this module found.', file=fd)
            print('\n Visit https://gitlab.com/ase/ase/tree/master/ase/io/{:}.py to see if you find the information needed in the source code.\n\n'.format(module_name), file=fd)


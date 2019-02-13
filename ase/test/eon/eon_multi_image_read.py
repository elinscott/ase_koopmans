
if False:
    import sys
    import os.path

    local_path = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir, os.path.pardir  ))

    print(local_path)
    b = sys.path
    sys.path = [local_path] + b


import ase.io

out_file = 'neb.con'

images = ase.io.read(out_file, format='eon')

print (len(images))

ase.io.write('neb_out.con', images, format='eon')

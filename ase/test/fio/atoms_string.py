from ase.build import bulk
from ase.io.string import atoms_to_string, atoms_to_bytes, parse_images


atoms = bulk('Ti')
images = [bulk('Au'), bulk('Ti'), bulk('NaCl', 'rocksalt', 17)]

txt = atoms_to_string(atoms, format='xyz')
print(txt)

txt1 = atoms_to_string(images, format='xyz')
print(txt1)

images1 = parse_images(txt1)
print(imgs)

print(images1)
print(images2)

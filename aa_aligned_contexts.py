#!/usr/bin/python

from __future__ import print_function

import __main__
__main__.pymol_argv = [ 'pymol', '-qc']

import sys
import os
#sys.path.insert(1, '/usr/local/Cellar/pymol/1.7.0.0/lib/python2.7/site-packages/')
import pymol
from pymol import cmd

pymol.finish_launching()

if len(sys.argv) != 6:
    print('Usage: %s pdb-file align-sele save-sele align-atoms frame-pdb')
    sys.exit(1)

struct_path = sys.argv[1]
struct_name = os.path.splitext(os.path.basename(struct_path))[0]
cmd.load(struct_path, struct_name)

cmd.split_states(struct_name, 1, 1, 'temp')
cmd.delete(struct_name)
# cmd.save('%s_debug.pse' % struct_name)
cmd.set_name('temp0001', struct_name)

frame_bname = 'frame'
sys
# frame_path = './aa_frames/ile_frame.pdb'
frame_path = sys.argv[5]
cmd.load(frame_path, frame_bname)

# target_sele = 'resn asp'
target_sele = sys.argv[2]
# dist = 6
#dist = sys.argv[3]
user_sele = sys.argv[3]

# 1) find residue number of all residues in the selection
# 2) align residue to the frame
# 3) create objects for every pair of interest
# 4) save object, delete

target_resis = set()
pymol_namespace = {"target_resis": target_resis}
cmd.iterate(target_sele, 'target_resis.add( (chain, resi, resn, alt) )', space=pymol_namespace)

target_resis = sorted(target_resis)

print('hits:')
print('chain, resi_icode, resn, altLoc')
print('\n'.join([str(x) for x in target_resis]))

# align_atoms = 'CA+CB+CG1'
align_atoms = sys.argv[4]
frame_align_sele = '/%s////%s' % (frame_bname, align_atoms)

secondary_sele_str = '%s and %s' % (struct_name, user_sele)

# serial number appended to names
serial = 0

for s in target_resis:

    print('\nloop iteration begin')

    chain = s[0]
    resi_icode = s[1]
    res_sele = '/'.join(['', struct_name, '', chain, resi_icode])

    target_align_sele = '%s/%s' % (res_sele, align_atoms)
    print('cmd used to select target resi =', target_align_sele)

    print('dumping begin debug')
    cmd.save('%s_debug_before_%04d.pse' % (struct_name, serial))

    # cmd.pair_fit(bb_sele, frame_sele)
    pair_fit_str = 'pair_fit %s, %s' % (target_align_sele, frame_align_sele)
    print(pair_fit_str)
    cmd.do(pair_fit_str)

    cmd.select('target', res_sele)
    # cmd.create('temp', secondary_sele_str)
    create_temp_str = 'create temp, %s' % secondary_sele_str
    print('temp object creation command =', create_temp_str)
    cmd.do(create_temp_str)
    print('done')

    print('dumping middle debug')
    cmd.save('%s_debug_middle_%04d.pse' % (struct_name, serial))

    #cmd.save('%s_%04d.pdb' % (struct_name, serial), 'temp')
    save_temp_str = 'save %s_%04d.pdb' % (struct_name, serial) + ', temp'
    print('save_temp_str =', save_temp_str)
    cmd.do(save_temp_str)
    print('done')

    cmd.save('%s_debug_after_%04d.pse' % (struct_name, serial))

    serial += 1
    #cmd.save(struct_name + '_' + chain + resi_icode + '_expand6.pse', 'all')
    cmd.delete('temp')







#def f(chain,resi_icode,resn,alt):
#    # if there's an insertion code it ends up in resi_icode
#    # ex: 180A
#    print "%s    %s    %s    %s" % (chain, resi_icode, resn, alt)
#
#namespace = {"f" : f, "sele_instances" : []}
#outname = "asp"
#sele_str = "resn asp and name CA"
#expand_dist = 3.5
#pymol.cmd.iterate(sele_str, "sele_instances.append( (chain, resi, resn, alt) )", space=namespace)
#
#i = 0;
#for sele_instances in namespace["sele_instances"]:
#    i += 1;
#    sele_and_context_sele = "byres ( ( c. " + sele_instances[0] + " and resi " + sele_instances[1] + ") expand " + str(expand_dist) + ")"
#    pymol.cmd.create("temp", sele_and_context_sele )
#
#    temp_frame_sele = "temp//" + sele_instances[0] + "/" + sele_instances[1] + "/N+CA+C"
#
#    pymol.cmd.pair_fit(temp_frame_sele, "aa_frame////N+CA+C")
#    pymol.cmd.save(structureName + "_" + outname + str(i) + "_expand3p5.pdb", "temp" )
#    pymol.cmd.delete("temp")

#pymol.cmd.quit()

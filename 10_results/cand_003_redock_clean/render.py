import sys
from pymol import cmd
R, name, cx_file, tag = sys.argv[1:5]
grp = "mao" if R.startswith("2") else "ttbk"
nat = {"2V5Z":"SAG","2Z5X":"HRM","7JXX":"VP7","7Q8Y":"9IV"}[R]
cmd.load(cx_file, "cx")
cmd.bg_color("white"); cmd.set("ray_opaque_background", 1); cmd.set("antialias", 2)
cmd.select("lig", "resn C03"); cmd.select("prot", "polymer.protein")
cmd.select("pocket", "byres (prot within 4.0 of lig)")
cmd.select("shell", "byres (prot within 11.0 of lig)")
cmd.hide("everything"); cmd.show("cartoon", "shell"); cmd.color("gray85", "shell"); cmd.set("cartoon_transparency", 0.6)
cmd.show("sticks", "pocket and not hydro"); cmd.set("stick_radius", 0.11, "pocket")
cmd.color("skyblue", "pocket and elem C"); cmd.util.cnc("pocket")
cmd.show("sticks", "lig and not (hydro and elem H and neighbor elem C)")
cmd.set("stick_radius", 0.24, "lig"); cmd.color("orange", "lig and elem C"); cmd.util.cnc("lig")
# crystal references, drawn as thin ghosts (they were REMOVED from the receptor that was docked)
cmd.load(f"../../03_receptors/{grp}/{R}/native_{nat}.pdb", "native")
cmd.show("lines", "native"); cmd.color("green", "native"); cmd.set("line_width", 2.5, "native")
if grp == "mao":
    cmd.load(f"../../03_receptors/{grp}/{R}/raw.pdb", "raw"); cmd.select("fad", "raw and chain A and resn FAD")
    cmd.hide("everything", "raw"); cmd.show("lines", "fad"); cmd.color("purple", "fad"); cmd.set("line_width", 2.5, "fad")
cmd.distance("hb", "lig and elem N+O and not hydro", "pocket and elem N+O and not hydro", 3.3, mode=0)
cmd.hide("labels", "hb"); cmd.color("black", "hb"); cmd.set("dash_gap", 0.3); cmd.set("dash_width", 3)
cmd.label("pocket and name CA and not resi 0", "resn+resi"); cmd.set("label_size", -0.7); cmd.set("label_color", "black")
cmd.set("label_position", (0, 0, 4)); cmd.set("depth_cue", 0); cmd.set("ray_shadows", 0)
cmd.orient("lig"); cmd.zoom("lig", 9)
cmd.png(f"img/{tag}.png", 1400, 1000, dpi=200, ray=1)

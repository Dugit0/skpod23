import os

base = 'base_prog'
gen = 'gen_prog'
N = [100, 200, 500, 1000]
for name in os.listdir(base):
    inp_name = os.path.join(base, name)
    with open(inp_name) as prog_file:
        lines = prog_file.readlines()
    for n in N:
        out_name = os.path.join(gen, f"{name.split('.')[0]}{n}.c")
        with open(out_name, "w") as f_out:
            for line in lines:
                if line.startswith("#define  N"):
                    f_out.write(f"#define  N  {n}\n")
                else:
                    f_out.write(line)

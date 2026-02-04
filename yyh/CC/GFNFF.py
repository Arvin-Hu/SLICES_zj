import os
import subprocess
os.environ["XTB_MOD_PATH"] = "/crystal/src/slices/xtb_noring_nooutput_nostdout_noCN"


def lmp_to_xyz(lmp_file, xyz_file):
    """将 LAMMPS 数据文件转换为 XYZ 格式（仅支持简单碳原子结构）"""
    atoms = []
    with open(lmp_file, 'r') as f:
        lines = f.readlines()
    atom_section = False
    for line in lines:
        if line.strip() == "Atoms":
            atom_section = True
            continue
        if atom_section:
            if line.strip() == "":
                break
            parts = line.strip().split()
            if len(parts) >= 6:
                # 格式: id type charge x y z
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                atoms.append((x, y, z))
    with open(xyz_file, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write("TBLG structure from LAMMPS data\n")
        for x, y, z in atoms:
            f.write(f"C {x:.6f} {y:.6f} {z:.6f}\n")

def optimize_with_xtb(xyz_file, out_file):
    """调用xtb进行GFN-FF优化"""
    cmd = f"/crystal/src/slices/xtb_noring_nooutput_nostdout_noCN {xyz_file} --gfnff"
    subprocess.run(cmd, shell=True, check=True)
    # xtb 优化后会生成 xtbopt.xyz
    if os.path.exists("xtbopt.xyz"):
        os.rename("xtbopt.xyz", out_file)

if __name__ == "__main__":
    lmp_file = "yyh/CC/TBLG_5.0deg.lmp"
    xyz_file = "TBLG_5.0deg.xyz"
    out_file = "TBLG_5.0deg_optimize.GFNFF"
    lmp_to_xyz(lmp_file, xyz_file)
    optimize_with_xtb(xyz_file, out_file)
    print(f"优化完成，结果已保存到 {out_file}")
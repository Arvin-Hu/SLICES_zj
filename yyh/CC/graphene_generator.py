import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import os

class BilayerTwistedGraphene:
    def __init__(self):
        # 石墨烯晶格常数
        self.a = 2.46  # 晶格常数 (Å)
        self.a0 = 1.42  # C-C键长 (Å)
        self.d_layer = 3.35  # 层间距 (Å)
        
        # 默认参数
        self.angle = 0.0  # 转角 (度)
        self.lx = 50.0  # 长方体x方向长度 (Å)
        self.ly = 50.0  # 长方体y方向长度 (Å)
        self.lz = 10.0  # 长方体z方向高度 (Å)
        
        # 原子数据
        self.atoms = []
        self.atom_types = []
        
    def generate_unit_cell(self, layer=1):
        """生成石墨烯单层原胞原子坐标"""
        # 石墨烯原胞基矢
        a1 = np.array([self.a * math.sqrt(3), 0, 0])
        a2 = np.array([self.a * math.sqrt(3) / 2, self.a * 3 / 2, 0])
        
        # 原胞内两个原子的相对位置
        atoms_pos = [
            np.array([0, 0, 0]),
            np.array([self.a * math.sqrt(3) / 2, self.a / 2, 0])
        ]
        
        return a1, a2, atoms_pos
    
    def generate_layer(self, layer_z, rotation_angle=0.0):
        """生成单层石墨烯"""
        a1, a2, basis = self.generate_unit_cell()
        
        # 计算需要的原胞数量以覆盖长方体区域
        n1 = int(math.ceil(self.lx / np.linalg.norm(a1))) + 2
        n2 = int(math.ceil(self.ly / np.linalg.norm(a2))) + 2
        
        layer_atoms = []
        
        # 生成超胞
        for i in range(-n1, n1 + 1):
            for j in range(-n2, n2 + 1):
                for atom in basis:
                    pos = i * a1 + j * a2 + atom
                    
                    # 应用旋转（绕z轴）
                    if rotation_angle != 0:
                        theta = math.radians(rotation_angle)
                        rot_matrix = np.array([
                            [math.cos(theta), -math.sin(theta), 0],
                            [math.sin(theta), math.cos(theta), 0],
                            [0, 0, 1]
                        ])
                        pos = np.dot(rot_matrix, pos)
                    
                    # 设置z坐标
                    pos[2] = layer_z
                    
                    # 检查是否在长方体区域内
                    if (abs(pos[0]) <= self.lx/2 and 
                        abs(pos[1]) <= self.ly/2 and 
                        abs(pos[2]) <= self.lz/2):
                        layer_atoms.append(pos)
        
        return layer_atoms
    
    def generate_structure(self, angle=None, lx=None, ly=None, lz=None):
        """生成双层转角石墨烯结构"""
        if angle is not None:
            self.angle = angle
        if lx is not None:
            self.lx = lx
        if ly is not None:
            self.ly = ly
        if lz is not None:
            self.lz = lz
        
        # 清空原子列表
        self.atoms = []
        self.atom_types = []
        
        # 生成底层（不旋转）
        layer1_atoms = self.generate_layer(-self.d_layer/2, 0)
        self.atoms.extend(layer1_atoms)
        self.atom_types.extend([1] * len(layer1_atoms))
        
        # 生成顶层（旋转指定角度）
        layer2_atoms = self.generate_layer(self.d_layer/2, self.angle)
        self.atoms.extend(layer2_atoms)
        self.atom_types.extend([2] * len(layer2_atoms))
        
        print(f"生成双层转角石墨烯结构:")
        print(f"  转角: {self.angle} 度")
        print(f"  长方体尺寸: {self.lx} × {self.ly} × {self.lz} Å")
        print(f"  原子总数: {len(self.atoms)}")
        print(f"  底层原子数: {len(layer1_atoms)}")
        print(f"  顶层原子数: {len(layer2_atoms)}")
        
        return np.array(self.atoms), self.atom_types
    
    def visualize(self, save_fig=False, filename="graphene_structure.png"):
        """可视化石墨烯结构"""
        atoms_array = np.array(self.atoms)
        if len(atoms_array) == 0:
            print("请先生成结构")
            return
        
        fig = plt.figure(figsize=(12, 5))
        
        # 2D顶视图
        ax1 = fig.add_subplot(121)
        layer1_mask = np.array(self.atom_types) == 1
        layer2_mask = np.array(self.atom_types) == 2
        
        ax1.scatter(atoms_array[layer1_mask, 0], atoms_array[layer1_mask, 1], 
                   c='blue', s=20, label='bottom', alpha=0.6)
        ax1.scatter(atoms_array[layer2_mask, 0], atoms_array[layer2_mask, 1], 
                   c='red', s=20, label='top', alpha=0.6)
        ax1.set_xlabel('X (Å)')
        ax1.set_ylabel('Y (Å)')
        ax1.set_title(f'TBLG (θ = {self.angle}°)')
        ax1.legend()
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        
        # 3D视图
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.scatter(atoms_array[layer1_mask, 0], atoms_array[layer1_mask, 1], 
                   atoms_array[layer1_mask, 2], c='blue', s=20, alpha=0.6)
        ax2.scatter(atoms_array[layer2_mask, 0], atoms_array[layer2_mask, 1], 
                   atoms_array[layer2_mask, 2], c='red', s=20, alpha=0.6)
        ax2.set_xlabel('X (Å)')
        ax2.set_ylabel('Y (Å)')
        ax2.set_zlabel('Z (Å)')
        ax2.set_title('3D view')
        
        plt.tight_layout()
        
        if save_fig:
            plt.savefig(filename, dpi=300)
            print(f"结构图已保存为 {filename}")
        
        plt.show()
    
    def save_lammps_data(self, filename="graphene.lmp"):
        """保存为LAMMPS数据文件，包含混合势函数配置信息"""
        if len(self.atoms) == 0:
            print("请先生成结构")
            return
        
        atoms_array = np.array(self.atoms)
        
        # 计算盒子边界
        x_min, x_max = atoms_array[:, 0].min(), atoms_array[:, 0].max()
        y_min, y_max = atoms_array[:, 1].min(), atoms_array[:, 1].max()
        z_min, z_max = atoms_array[:, 2].min(), atoms_array[:, 2].max()
        
        # 添加一些真空区域
        padding = 10.0  # 增加padding确保边界足够
        x_min -= padding
        x_max += padding
        y_min -= padding
        y_max += padding
        z_min -= padding
        z_max += padding
        
        with open(filename, 'w') as f:
            # ==================== 文件头信息 ====================
            f.write("#" + "=" * 60 + "\n")
            f.write(f"# 双层转角石墨烯 LAMMPS 数据文件\n")
            f.write(f"# 转角角度: {self.angle}°\n")
            f.write(f"# 结构尺寸: {self.lx:.1f} × {self.ly:.1f} × {self.lz:.1f} Å\n")
            f.write(f"# 层间距: {self.d_layer:.2f} Å\n")
            f.write(f"# 总原子数: {len(self.atoms)}\n")
            f.write(f"# 底层原子数: {np.sum(np.array(self.atom_types) == 1)}\n")
            f.write(f"# 顶层原子数: {np.sum(np.array(self.atom_types) == 2)}\n")
            f.write("#" + "=" * 60 + "\n\n")
            
            # ==================== 基本数据部分 ====================
            f.write(f"{len(self.atoms)} atoms\n")
            f.write("2 atom types\n\n")
            
            f.write(f"{x_min:.6f} {x_max:.6f} xlo xhi\n")
            f.write(f"{y_min:.6f} {y_max:.6f} ylo yhi\n")
            f.write(f"{z_min:.6f} {z_max:.6f} zlo zhi\n\n")
            
            # ==================== 原子质量 ====================
            f.write("Masses\n\n")
            f.write("# 原子类型 质量(amu) 注释\n")
            f.write("1 12.0107  # 底层碳原子\n")
            f.write("2 12.0107  # 顶层碳原子\n\n")
            
            # ==================== 原子坐标 ====================
            f.write("Atoms # atomic\n\n")
            f.write("# id atom_type charge x y z\n")
            for i, (atom, atom_type) in enumerate(zip(self.atoms, self.atom_types), 1):
                f.write(f"{i} {atom_type} 0 {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n")
            
            # ==================== 邻接信息（可选） ====================
            f.write("\n# 注：邻接信息通常由LAMMPS在运行时自动生成\n")
            f.write("# 如果需要手动指定，可以在此添加Bonds、Angles、Dihedrals等部分\n")
        
        print(f"LAMMPS数据文件已保存为 {filename}")
        
        # ==================== 自动生成配套的输入脚本 ====================
        self._generate_lammps_input_script(filename)
        
    def _generate_lammps_input_script(self, data_filename):
        """生成与数据文件配套的LAMMPS输入脚本"""
        script_name = data_filename.replace(".lmp", "_input.in")
        
        with open(script_name, 'w') as f:
            f.write("#" + "=" * 70 + "\n")
            f.write("# 双层转角石墨烯 LAMMPS 输入脚本\n")
            f.write("# 自动生成 - 配合数据文件: " + data_filename + "\n")
            f.write("#" + "=" * 70 + "\n\n")
            
            f.write("# ============ 基本设置 ============\n")
            f.write("units metal\n")
            f.write("atom_style atomic\n")
            f.write("boundary p p p\n")
            f.write("atom_modify map array\n\n")
            
            f.write("# ============ 读取数据文件 ============\n")
            f.write(f"read_data {data_filename}\n\n")
            
            f.write("# ============ 势函数设置 ============\n")
            f.write("# 使用 hybrid/overlay 组合两种势函数\n")
            f.write("# AIREBO 用于层内共价键，LJ 用于层间范德华相互作用\n")
            f.write("pair_style hybrid/overlay airebo 3.0 lj/cut 15.0\n\n")
            
            f.write("# 设置AIREBO势函数参数：用于层内相互作用\n")
            f.write("pair_coeff 1 1 airebo CH.airebo C C  # 底层-底层\n")
            f.write("pair_coeff 2 2 airebo CH.airebo C C  # 顶层-顶层\n\n")
            
            f.write("# 设置LJ势函数参数：用于层间相互作用\n")
            f.write("# 参数来源: Kolmogorov, A. N. & Crespi, V. H. Phys. Rev. B 71, 235415 (2005)\n")
            f.write("# epsilon = 2.39 meV = 0.00239 eV, sigma = 3.41 Å\n")
            f.write("pair_coeff 1 2 lj/cut 0.00239 3.41\n\n")
            
            f.write("# 注：对于同一类型原子之间的非键相互作用，AIREBO已包含\n")
            f.write("# 因此不需要再为1-1或2-2对设置LJ势\n\n")
            
            f.write("# ============ 邻居列表设置 ============\n")
            f.write("neighbor 2.0 bin\n")
            f.write("neigh_modify every 10 delay 0 check yes\n\n")
            
            f.write("# ============ 能量最小化 ============\n")
            f.write("thermo 100\n")
            f.write("thermo_style custom step temp pe ke etotal press vol\n")
            f.write("minimize 1.0e-10 1.0e-10 10000 10000\n\n")
            
            f.write("# ============ 保存最小化后的结构 ============\n")
            f.write("write_data minimized.lmp\n")
            f.write("write_restart minimized.restart\n\n")
            
            f.write("# ============ 可选：可视化输出 ============\n")
            f.write("dump 1 all atom 1000 relaxed.xyz\n")
            f.write("dump_modify 1 element C\n\n")
            
            f.write("# ============ 可选：NVT动力学模拟 ============\n")
            f.write("# 取消注释以下行以运行动力学模拟\n")
            f.write("# reset_timestep 0\n")
            f.write("# velocity all create 300.0 12345\n")
            f.write("# fix 1 all nvt temp 300.0 300.0 0.1\n")
            f.write("# timestep 0.001\n")
            f.write("# run 10000\n\n")
            
            f.write("# ============ 可选：层间相互作用分析 ============\n")
            f.write("# 计算层间相互作用能\n")
            f.write("# group layer1 type 1\n")
            f.write("# group layer2 type 2\n")
            f.write("# compute inter layer2 group/group layer1\n")
            f.write("# variable inter_energy equal c_inter\n")
            f.write("# run 0\n")
            f.write("# print \"层间结合能 = ${inter_energy} eV\"\n\n")
            
            f.write("# ============ 脚本结束 ============\n")
            f.write("print \"模拟完成！\"\n")
        
        print(f"配套的LAMMPS输入脚本已生成: {script_name}")
        print("=" * 60)
        print("重要提醒：")
        print("1. 请确保CH.airebo势函数文件在运行目录中")
        print("2. 可以从https://github.com/lammps/lammps/tree/master/potentials下载势函数文件")
        print("3. 首次运行前建议调整参数以适应您的具体需求")
        print("=" * 60)
        
        def save_xyz(self, filename="graphene.xyz"):
            """保存为XYZ格式文件"""
            if len(self.atoms) == 0:
                print("请先生成结构")
                return
            
            with open(filename, 'w') as f:
                f.write(f"{len(self.atoms)}\n")
                f.write(f"双层转角石墨烯，转角={self.angle}°，尺寸={self.lx}x{self.ly}x{self.lz}Å\n")
                
                for atom, atom_type in zip(self.atoms, self.atom_types):
                    element = 'C'
                    f.write(f"{element} {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n")
            
            print(f"XYZ文件已保存为 {filename}")


def main():
    """主函数：演示如何使用"""
    print("=" * 60)
    print("双层转角石墨烯结构生成器")
    print("=" * 60)
    
    # 创建生成器实例
    generator = BilayerTwistedGraphene()
    
    # 示例1：小角度转角石墨烯（魔角附近）
    print("\n示例1: 生成小角度转角石墨烯 (θ = 1.05°)")
    generator.generate_structure(angle=1.05, lx=80, ly=80, lz=15)
    generator.visualize(save_fig=True, filename="twisted_graphene_1.05deg.png")
    generator.save_lammps_data("twisted_graphene_1.05deg.lmp")
    generator.save_xyz("twisted_graphene_1.05deg.xyz")
    
    # 示例2：较大角度转角
    print("\n\n示例2: 生成较大角度转角石墨烯 (θ = 10°)")
    generator.generate_structure(angle=10.0, lx=60, ly=60, lz=12)
    generator.visualize(save_fig=True, filename="twisted_graphene_10deg.png")
    generator.save_lammps_data("twisted_graphene_10deg.lmp")
    
    # 示例3：零角度（AA堆叠）
    print("\n\n示例3: 生成零转角石墨烯 (AA堆叠)")
    generator.generate_structure(angle=0.0, lx=40, ly=40, lz=10)
    generator.visualize(save_fig=True, filename="graphene_AA.png")
    generator.save_lammps_data("graphene_AA.lmp")
    
    # 示例4：AB堆叠（60度转角）
    print("\n\n示例4: 生成AB堆叠石墨烯 (θ = 60°)")
    generator.generate_structure(angle=60.0, lx=50, ly=50, lz=10)
    generator.visualize(save_fig=True, filename="graphene_AB.png")
    generator.save_lammps_data("graphene_AB.lmp")
    
    print("\n" + "=" * 60)
    print("所有示例生成完成！")
    print("=" * 60)


def interactive_mode():
    """交互模式：用户自定义参数"""
    print("=" * 60)
    print("双层转角石墨烯结构生成器 - 交互模式")
    print("=" * 60)
    
    generator = BilayerTwistedGraphene()
    
    # 获取用户输入
    try:
        angle = float(input("请输入转角角度 (度): "))
        lx = float(input("请输入长方体x方向长度 (Å): "))
        ly = float(input("请输入长方体y方向长度 (Å): "))
        lz = float(input("请输入长方体z方向高度 (Å): "))
        
        # 生成结构
        generator.generate_structure(angle=angle, lx=lx, ly=ly, lz=lz)
        
        # 可视化
        generator.visualize(save_fig=True)
        
        # 保存文件
        save_lmp = input("是否保存为LAMMPS数据文件? (y/n): ").lower()
        if save_lmp == 'y':
            filename = input("请输入文件名 (默认: graphene.lmp): ")
            if not filename:
                filename = "graphene.lmp"
            generator.save_lammps_data(filename)
        
        save_xyz = input("是否保存为XYZ文件? (y/n): ").lower()
        if save_xyz == 'y':
            filename = input("请输入文件名 (默认: graphene.xyz): ")
            if not filename:
                filename = "graphene.xyz"
            generator.save_xyz(filename)
            
    except ValueError:
        print("输入错误！请确保输入的是数字。")


if __name__ == "__main__":
    # 运行主演示
    main()
    
    # 或者运行交互模式
    # interactive_mode()
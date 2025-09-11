from typing import List, Literal, Optional, Union, Dict
from pydantic import BaseModel, Field, model_validator, confloat, conint, ValidationError

# ===== Core Sections =====

class Calculation(BaseModel):
    """&calculation / 计算理论与模式"""
    theory: Literal['dft', 'dft_md', 'tddft_response', 'tddft_pulse'] = 'dft'

class Control(BaseModel):
    """&control / 表示信息，不参与计算"""
    sysname: str = "System"
    comment: Optional[str] = None

class Unit(BaseModel):
    """&unit / 单位"""
    unit_system: Literal['A_eV_fs', 'a.u.'] = 'A_eV_fs'

class System(BaseModel):
    """&system / 体系规模，周期性，电子数等"""
    nelem: conint(ge=1)
    natom: conint(ge=1)
    nelec: conint(ge=1)
    nstate: conint(ge=1)
    yn_periodic: Literal['y', 'n'] = 'n'

    @model_validator(mode='after')
    def check_state(self):
        if self.nstate < (self.nelec + 1) // 2:
            raise ValueError('nstate must be >= ceil(nelec/2) as a basic guard')
        return self

class Pseudopotential(BaseModel): 
    """单个赝势元素的数据结构"""
    file_pseudo: str
    izatom: conint(ge=1)
    lloc_ps: conint(ge=0)
    index: conint(ge=1) # 添加 index 字段用于跨模型验证

class Pseudo(BaseModel): 
    """&pseudo / 赝势（按元素序 i=1..nelem）"""
    pseudopots: List[Pseudopotential] = Field(min_length=1)

class Functional(BaseModel):
    """&functional / 交换相关函数"""
    xc: Literal['PZ', 'PZM', 'TBmBJ'] = 'PZ'

class RGrid(BaseModel):
    """&rgrid / 实空间网格（孤立分子常用）"""
    dl: List[confloat(gt=0)] = Field(min_length=3, max_length=3, default=[0.25, 0.25, 0.25])
    num_rgrid: List[conint(ge=8)] = Field(min_length=3, max_length=3, default=[64, 64, 64])

class TGrid(BaseModel):
    """&tgrid / 动力学时间网格"""
    dt: confloat(gt=0) = 1.25e-3
    nt: conint(ge=1) = 5000

class Scf(BaseModel):
    """&scf / 电子结构计算相关参数"""
    nscf: conint(ge=1) = 300
    threshold: confloat(gt=0) = 1e-17

# ===== Emfield: 判别联合，按 shape 决定必填项 =====
class EmfieldImpulse(BaseModel):
    ae_shape1: Literal['impulse'] = 'impulse'
    epdir_re1: List[float] = Field(min_length=3, max_length=3, default=[1.0, 0.0, 0.0])

    @model_validator(mode='after')
    def _norm_check(self):
        import math
        n = math.sqrt(sum(d*d for d in self.epdir_re1))
        if n == 0:
            raise ValueError("epdir_re1 cannot be zero vector")
        return self

class EmfieldAcos(BaseModel):
    ae_shape1: Literal['Acos2', 'Acos3', 'Acos4']
    epdir_re1: List[float] = Field(min_length=3, max_length=3, default=[1.0, 0.0, 0.0])
    I_wcm2_1: confloat(ge=0) = 0.0
    tw1: confloat(gt=0) = 1.0
    omega1: confloat(ge=0) = 0.0

    @model_validator(mode='after')
    def _norm_check(self):
        import math
        n = math.sqrt(sum(d*d for d in self.epdir_re1))
        if n == 0:
            raise ValueError("epdir_re1 cannot be zero vector")
        return self

Emfield = Union[EmfieldImpulse, EmfieldAcos]

# ===== Analysis: 判别联合，按 dft/tddft 决定必填项 =====
class Analysisdft(BaseModel):
    """&analysis / 计算结果分析"""
    yn_out_psi:  Literal['y','n'] = 'y'
    yn_out_dns:  Literal['y','n'] = 'y'
    yn_out_dos:  Literal['y','n'] = 'y'
    yn_out_pdos: Literal['y','n'] = 'y'
    yn_out_elf:  Literal['y','n'] = 'y'

class Analysistddft(BaseModel):
    de: confloat(ge=0) = 0.01
    nenergy: conint(ge=1) = 1000

Analysis = Union[Analysisdft, Analysistddft]

class Atomic_Core(BaseModel):
    """&atomic_core / 原子核与坐标（最后一列是元素序 index ∈ [1..nelem]）"""
    class Atom(BaseModel):
        symbol: str
        x: float
        y: float
        z: float
        index: conint(ge=1)
    atoms: List[Atom] = Field(min_length=1)


# ===== 汇总：CommonCore =====

class CommonCore(BaseModel):
    calculation: Calculation = Calculation()
    control: Control = Control()
    unit: Unit = Unit()
    system: System
    pseudo: Pseudo
    functional: Functional = Functional()
    rgrid: Optional[RGrid] = None
    tgrid: Optional[TGrid] = None
    scf: Optional[Scf] = None
    emfield: Optional[Emfield] = None
    analysis: Analysis = Field(default_factory = Analysisdft)
    atomic_core: Atomic_Core

    @model_validator(mode='after')
    def check_consistency(self):
        # 1) natom 必须等于原子数
        if self.system.natom != len(self.atomic_core.atoms):
            raise ValueError(f"system.natom ({self.system.natom}) != number of atoms ({len(self.atomic_core.atoms)})")
            
        # 2) 赝势的元素数必须与 system.nelem 相同
        nelem = self.system.nelem
        if len(self.pseudo.pseudopots) != nelem:
            raise ValueError(f"pseudo.pseudopots 列表长度 ({len(self.pseudo.pseudopots)}) 必须与 system.nelem ({nelem}) 相同")

        # 3) 原子 index 必须在 [1..nelem] 范围内，并且与赝势索引一致
        atom_indices = set(a.index for a in self.atomic_core.atoms)
        pseudo_indices = set(p.index for p in self.pseudo.pseudopots)
        
        if not atom_indices.issubset(pseudo_indices):
            # 原子索引必须是赝势索引的子集，即每个原子都必须有对应的赝势定义
            missing_indices = sorted(list(atom_indices - pseudo_indices))
            raise ValueError(f"以下原子索引在赝势列表中找不到：{missing_indices}")
            
        if not pseudo_indices.issubset(atom_indices):
            # 每个赝势定义也必须至少被一个原子引用
            extra_indices = sorted(list(pseudo_indices - atom_indices))
            raise ValueError(f"以下赝势索引在原子列表中没有被引用：{extra_indices}")

        # 4) 检查原子索引的连续性和最大值
        if len(pseudo_indices) != nelem or max(pseudo_indices) != nelem:
            # 这里的逻辑可以根据你的具体需求调整。
            # 如果要求所有 1..nelem 的索引都必须被使用，可以这样检查。
            expected_indices = set(range(1, nelem + 1))
            if pseudo_indices != expected_indices:
                raise ValueError(f"赝势索引必须是连续的 [1..{nelem}]")
        
        return self
import math
from typing import List, Tuple, Any, Dict
PACK_GROUND = "ground_state"
PACK_POLAR = "polarizability"
PACK_ELEC_DYN = "electron_dynamics"

REQUIRED_BY_PROFILE = {
    "isolated" : [
        "core.calculation.theory",
        "core.control.sysname",
        "core.unit.unit_system",
        "core.system.nelem",
        "core.system.natom",
        "core.system.nelec",
        "core.system.nstate",
        "core.system.yn_periodic",
        "core.pseudo.pseudopots",
        "core.functional.xc",
        "core.atomic_core.atoms",
        "core.rgrid.dl",
        "core.rgrid.num_grid",
    ],
    "periodic" : [
       
    ],
    "mutisclale" : [
        
    ],
    "fdtd": [
        
    ]
   
}    

REQUIRED_BY_PACK = {
    PACK_GROUND : [
        "core.scf.nscf",
        "core.scf.threshold",
    ],
    PACK_POLAR : [
        "core.tgrid.dt","core.tgrid.nt",
        "core.emfield.as=shape1", "core.emfield.epdir_re1",
        "core.analysis.de",
        "core.analysis.nenergy"
    ],
    PACK_ELEC_DYN: [
        "core.tgrid.dt","core.tgrid.nt",
        "core.emfield.ae_shape1","core.emfield.I_wcm2_1","core.emfield.tw1","core.emfield.omega1","core.emfield.epdir_re1",
        "core.analysis.de",
        "core.analysis.nenergy"
    ]    
}

DEFAULTS = {
    "core.functional.xc": "PZ",
    "core.rgrid.dl": [0.25,0.25,0.25],
    "core.rgrid.num_rgrid": [64,64,64],
    "core.tgrid.dt": 1.25e-3,
    "core.tgrid.nt": 5000,
    # ground state defaults:
    f"{PACK_GROUND}core.scf.nscf": 300,
    f"{PACK_GROUND}core.scf.threshold": 1.0e-9,
    f"{PACK_GROUND}core.analysis.yn_out_psi": "y",
    f"{PACK_GROUND}core.analysis.yn_out_dns": "y",
    f"{PACK_GROUND}core.analysis.yn_out_dos": "y",
    f"{PACK_GROUND}core.analysis.yn_out_pdos": "y",
    f"{PACK_GROUND}core.analysis.yn_out_elf": "y",
    # Polar defaults:
    f"{PACK_POLAR}.core.emfield.ae_shape1": "impulse",
    f"{PACK_POLAR}.core.emfield.epdir_re1": [0,0,1],
    f"{PACK_POLAR}.core.analysis.de": 0.01,
    f"{PACK_POLAR}.core.analysis.nenergy": 5000,
    # elecdyn defaults:
    f"{PACK_ELEC_DYN}.core.emfield.ae_shape1": "Acos2",
    f"{PACK_ELEC_DYN}.core.emfield.I_wcm2_1": 5.0e13,
    f"{PACK_ELEC_DYN}.core.emfield.tw1": 6.0,
    f"{PACK_ELEC_DYN}.core.emfield.omega1": 3.0,
    f"{PACK_ELEC_DYN}.core.emfield.epdir_re1": [0,0,1],
    f"{PACK_ELEC_DYN}.core.analysis.de": 0.01,
    f"{PACK_ELEC_DYN}.core.analysis.nenergy": 5000
}

def plan(spec: dict) -> dict:
    #只支持isolated molecule
    profile = spec.get("profile","isolated")
    if profile not in REQUIRED_BY_PROFILE:
        return {"success": False, "enabled_blocks": [], "warings": [],
                "message": "Sorry, only 'isolated molecule' simulation is supported now. Please wait for follwing updates."}
    
    #归一化packs，默认情况为计算基态
    pack = [p.get("name") for p in spec.get("packs", [])]
    if not packs and spec.get("core",{}). get("calculation",{}).get("theory","dft") == "dft":
        pack = [PACK_GROUND]

    #启用blocks
    enabled = ["calculation","control","unit","system","pseudo","functional","atomic_core","rgrid"]
    if PACK_GROUND in pack:
        enabled.append("scf")
        enabled.append("analysis_dft")
    if PACK_POLAR in pack:
        enabled += ["tgrid_polar","emfield_impulse","analysis_tddft"]
    if PACK_ELEC_DYN in pack:
        enabled += ["tgrid_edyn","emfield_acos","analysis_tddft"]

    seen = set(); ebaled_blocks = [b for b in enabled if not (b in seen or seen.add(b))]
        
    #缺少参数与默认建议
    missing = BASE_REQUIRED_ISOLATED.copy()
    for pk in packs:
        missing += REQUIRED_BY_PACK.get(pk, [])
    # 计算缺参
    def _get(d, path):
        cur=d
        for k in path.split("."):
            if not isinstance(cur, dict) or k not in cur: return None
            cur=cur[k]
        return cur
    missing_params = [path for path in missing if _get(spec, path) is None]

    # 建议默认：通用 + 仅选中 packs 的默认
    suggested = {}
    def _set(d, path, val):
        cur=d; ks=path.split(".")
        for k in ks[:-1]:
            cur.setdefault(k, {})
            cur=cur[k]
        cur[ks[-1]]=val

    # 通用默认
    for k,v in DEFAULTS.items():
        if any(k.startswith(p+".") for p in (PACK_GROUND,PACK_POLAR,PACK_EDYN)):
            continue
        if _get(spec, k) is None:
            _set(suggested, k, v)
    # 选中的 packs 默认
    for pk in packs:
        prefix = pk + "."
        for k,v in DEFAULTS.items():
            if k.startswith(prefix):
                key = k[len(prefix):]
                if _get(spec, key) is None:
                    _set(suggested, key, v)

    # 合并形成 proposed_payload（不覆盖用户已有）
    proposed = {"profile": profile, "packs": spec.get("packs", []), "core": spec.get("core", {})}
    for path in missing_params:
        val = _get(suggested, path)
        if val is not None and _get(proposed, path) is None:
            _set(proposed, path, val)

    checklist = [
        f"场景: {profile}", f"功能: {[p for p in packs] or '无'}", "—— 缺少的关键参数 ——"
    ] + [
        (f"[缺] {m} → 建议: { _get(suggested,m) }" if _get(suggested,m) is not None else f"[缺] {m}")
        for m in missing_params
    ] + ["是否使用以上建议默认来补齐？"]

    return {
        "success": True,
        "enabled_blocks": enabled_blocks,
        "missing_params": missing_params,
        "suggested_defaults": suggested,
        "proposed_payload": proposed,
        "checklist": checklist,
        "warnings": [],
    }
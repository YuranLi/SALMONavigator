class InputRenderer:
    def __init__(self):
        self.output = []
    def _add_line(self, line: str = ""):
        self.output.append(line)
    def _render_common_header(self, title: str):
        self._add_line("!########################################################################################!")
        self._add_line(f"! Object : {title} of {self.model.control.sysname} molecule")
        self._add_line("!----------------------------------------------------------------------------------------!")
        self._add_line("! * Conversion from unit_system = 'a.u.' to 'A_eV_fs':                                   !")
        self._add_line("!   Length: 1 [a.u.] = 0.52917721067    [Angstrom]                                       !")
        self._add_line("!   Energy: 1 [a.u.] = 27.21138505      [eV]                                             !")
        self._add_line("!   Time  : 1 [a.u.] = 0.02418884326505 [fs]                                             !")
        self._add_line("!########################################################################################!")
    def _render_core_sections(self):
        self._add_line(f"&calculation")
        self._add_line(f"  theory = {repr(self.model.calculation.theory)}")
        self._add_line(f"/")
        self._add_line("")
        self._add_line("&control")
        self._add_line(f"  sysname = {repr(self.model.control.sysname)}")
        self._add_line(f"/")
        self._add_line("")
        self._add_line("&units")
        self._add_line(f"  unit_system = {repr(self.model.unit.unit_system)}")
        self._add_line("/")
        self._add_line("")
        self._add_line("&system")
        self._add_line(f"  yn_periodic = {repr(self.model.system.yn_periodic)}")
        self._add_line(f"  nelem = {self.model.system.nelem}")
        self._add_line(f"  natom = {self.model.system.natom}")
        self._add_line(f"  nelec = {self.model.system.nelec}")
        self._add_line(f"  nstate = {self.model.system.nstate}")
        self._add_line("/")
        self._add_line("")
        self._add_line("&pseudo")
        for pot in self.model.pseudo.pseudopots:
            self._add_line(f"  file_pseudo({pot.index}) = {repr(pot.file_pseudo)}")
            self._add_line(f"  izatom({pot.index}) = {pot.izatom}")
            self._add_line(f"  lloc_ps({pot.index}) = {pot.lloc_ps}")
        self._add_line("/")
        self._add_line("")
        self._add_line("&functional")
        self._add_line(f"  xc = {repr(self.model.functional.xc)}")
        self._add_line("/")
        self._add_line("")
        if self.model.rgrid:
            self._add_line("&rgrid")
            dl_str = ', '.join([f'{d:.2f}d0' for d in self.model.rgrid.dl])
            num_str = ', '.join([str(n) for n in self.model.rgrid.num_rgrid])
            self._add_line(f"  dl(1:3) = {dl_str}")
            self._add_line(f"  num_rgrid(1:3) = {num_str}")
            self._add_line("/")
            self._add_line("")
    def _render_atomic_coor(self):
        self._add_line("&atomic_coor")
        for atom in self.model.atomic_core.atoms:
            self._add_line(f"  {repr(atom.symbol):<5} {atom.x:>10.6f} {atom.y:>10.6f} {atom.z:>10.6f}  {atom.index}")
        self._add_line("/")
        self._add_line("")
    def render_ground_state(self, model: CommonCore) -> str:
        self.output = []
        self.model = model
        self._render_common_header("Ground state")
        self._add_line("")
        self._render_core_sections()
        if self.model.scf:
            self._add_line("&scf")
            self._add_line(f"  nscf = {self.model.scf.nscf}")
            self._add_line(f"  threshold = {self.model.scf.threshold:.1e}")
            self._add_line("/")
            self._add_line("")
        if isinstance(self.model.analysis, Analysisdft):
            self._add_line("&analysis")
            self._add_line(f"  yn_out_psi = {repr(self.model.analysis.yn_out_psi)}")
            self._add_line(f"  yn_out_dns = {repr(self.model.analysis.yn_out_dns)}")
            self._add_line(f"  yn_out_dos = {repr(self.model.analysis.yn_out_dos)}")
            self._add_line(f"  yn_out_pdos = {repr(self.model.analysis.yn_out_pdos)}")
            self._add_line(f"  yn_out_elf = {repr(self.model.analysis.yn_out_elf)}")
            self._add_line("/")
            self._add_line("")
        self._render_atomic_coor()
        return "\n".join(self.output)
    def render_polarizability(self, model: CommonCore) -> str:
        self.output = []
        self.model = model
        self._render_common_header("Polarizability and photoabsorption")
        self._add_line("")
        self._render_core_sections()
        if self.model.tgrid:
            self._add_line("&tgrid")
            self._add_line(f"  dt = {self.model.tgrid.dt}")
            self._add_line(f"  nt = {self.model.tgrid.nt}")
            self._add_line("/")
            self._add_line("")
        if self.model.emfield and isinstance(self.model.emfield, EmfieldImpulse):
            self._add_line("&emfield")
            self._add_line(f"  ae_shape1 = {repr(self.model.emfield.ae_shape1)}")
            epdir_str = ', '.join([str(x) for x in self.model.emfield.epdir_re1])
            self._add_line(f"  epdir_re1(1:3) = {epdir_str}")
            self._add_line("/")
            self._add_line("")
        if isinstance(self.model.analysis, Analysistddft):
            self._add_line("&analysis")
            self._add_line(f"  de = {self.model.analysis.de}")
            self._add_line(f"  nenergy = {self.model.analysis.nenergy}")
            self._add_line("/")
            self._add_line("")
        self._render_atomic_coor()
        return "\n".join(self.output)
    def render_electron_dynamics(self, model: CommonCore) -> str:
        self.output = []
        self.model = model
        self._render_common_header("Electron dynamics")
        self._add_line("")
        self._render_core_sections()
        if self.model.tgrid:
            self._add_line("&tgrid")
            self._add_line(f"  dt = {self.model.tgrid.dt}")
            self._add_line(f"  nt = {self.model.tgrid.nt}")
            self._add_line("/")
            self._add_line("")
        if self.model.emfield and isinstance(self.model.emfield, EmfieldAcos):
            self._add_line("&emfield")
            self._add_line(f"  ae_shape1 = {repr(self.model.emfield.ae_shape1)}")
            self._add_line(f"  I_wcm2_1 = {self.model.emfield.I_wcm2_1:.2e}")
            self._add_line(f"  tw1 = {self.model.emfield.tw1}")
            self._add_line(f"  omega1 = {self.model.emfield.omega1}")
            epdir_str = ', '.join([str(x) for x in self.model.emfield.epdir_re1])
            self._add_line(f"  epdir_re1(1:3) = {epdir_str}")
            self._add_line("/")
            self._add_line("")
        if isinstance(self.model.analysis, Analysistddft):
            self._add_line("&analysis")
            self._add_line(f"  de = {self.model.analysis.de}")
            self._add_line(f"  nenergy = {self.model.analysis.nenergy}")
            self._add_line("/")
            self._add_line("")
        self._render_atomic_coor()
        return "\n".join(self.output)


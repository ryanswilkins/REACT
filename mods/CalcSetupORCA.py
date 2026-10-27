#
# Copyright (C) 2022 Geir Villy Isaksen and Bente Sirin Barge
#
# This file is part of chemREACT.
#
# chemREACT is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# chemREACT is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along
# with chemREACT. If not, see <https://www.gnu.org/licenses/>.
#

import os
from os import path, mkdir, remove
import glob
import shutil
from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSlot, QTimer, Qt
from UIs.SetupWindowOrca import Ui_SetupWindow
from mods.ScanBond import AtomBond
from mods.common_functions import atom_distance, random_color, write_file
import copy

# Short, read-only example content shown in the %Block preview panel when a
# block type is selected (see update_block_preview) - the inner content only
# (what goes between "%blockname" and "end"), for reference, not inserted
# automatically. Populated with real ORCA manual syntax where verified;
# deliberately left out of this dict (falling back to the generic "see the
# ORCA manual" message in update_block_preview) for anything not confirmed
# against the manual.
# ("mrcc" and "casresp" are intentionally omitted - verified keyword names
# exist for casresp, but not a confirmed complete/mandatory combination;
# mrcc's syntax couldn't be confirmed at all.)
BLOCK_KEYWORD_PREVIEWS = {
    "autoci": "citype FIC-MRCI\n"
    "# requires a preceding %casscf block (nel/norb/mult/nroots)",
    "basis": 'Basis "def2-TZVP"       # orbital basis\n'
    'AuxJ  "def2/J"          # RI-J auxiliary basis\n'
    'AuxC  "def2-TZVP/C"     # auxiliary basis for correlated methods',
    "casscf": "nel  10   # active-space electrons\n"
    "norb  8   # active-space orbitals\n"
    "mult  1   # multiplicity",
    "cim": "CIMTHRESH 0.0005   # default is 0.001",
    "cis": "NRoots 10\nTriplets true",
    "coords": "CTyp xyz\nCharge 0\nMult 1\ncoords\n"
    "  O   0.000   0.000   0.000\n"
    "  H   0.000   0.000   0.960\n"
    "  H   0.927   0.000  -0.240\n"
    "end",
    "compound": "New_Step\n! B3LYP def2-SVP Opt\n* xyzfile 0 1 start.xyz\n"
    "Step_End\nEndRun",
    "cosmors": 'solvent "Water"   # picks a built-in COSMO-RS solvent',
    "cpcm": "surfacetype gepol_ses_gaussian   # vdw_gaussian (default), gepol_ses, "
    "gepol_sas, gepol_ses_gaussian\n"
    "radius[8]   1.3                  # custom atomic radius (A) for element "
    "8 (oxygen)\n"
    "# epsilon/refrac/smd/draco/cds_cpcm are already handled by the "
    "Continuum solvation controls above",
    "elprop": "Dipole true\nPolar true",
    "eprnmr": "gtensor true\nNuclei = all H { aiso, adip }",
    "esd": 'GSHESSIAN "mymol.hess"\nHESSFLAG VG',
    "freq": "Temp      298.15   # temperature (K) for thermochemistry\n"
    "Increment 0.005    # numerical displacement step (bohr), if NumFreq",
    "geom": "Constraints\n"
    "  { C 0 C }      # freeze atom 0 at its current Cartesian position\n"
    "end",
    "ice": "nel     8       # active electrons\n"
    "norb    23      # active orbitals\n"
    "nroots  1\n"
    "mult    1\n"
    "tgen    1e-04   # generator threshold",
    "irc": "MaxIter    20\nDirection  both\nInitHess   read",
    "loc": "LocMet PM   # Pipek-Mezey localization",
    "mcrpa": "NRoots 8   # number of excited states\n"
    "# requires a preceding ground-state %casscf block (nroots 1)",
    "md": "Timestep    0.5_fs\n"
    "Initvel     300_K\n"
    "Thermostat  NHC 300_K Timecon 10.0_fs\n"
    "Run         100",
    "mdci": "citype  CCSD\nTriples 1          # add perturbative (T) correction",
    "method": "RI on   # use the RI-J approximation",
    "mp2": "RI       true\n"
    "EMin     -1.5     # orbital energy cutoff defining the frozen core (Eh)\n"
    "Density  relaxed  # relaxed density, needed for properties/gradients",
    "mrci": "CIType       SORCI\n"
    "NewBlock 1 0        # multiplicity 1, irrep 0\n"
    "  NRoots 1\n"
    "  Excitations cisd\n"
    "  Refs CAS(2,2) end\n"
    "end",
    "nbo": 'NBOKEYLIST = "$NBO NBO NPA AONBO=C ARCHIVE $END"   '
    "# requires NBOEXE env var",
    "neb": 'NEB_END_XYZFILE  "<change to product>.xyz"   # mandatory - end-point structure\n'
    'NEB_TS_XYZFILE "<change to TS guess>.xyz" #optional - if you have a TS guess\n'
    "NImages  8",
    "numgrad": "CentralDiff true   # two-sided differentiation\n"
    "DX          0.005   # displacement increment (bohr)",
    "output": "PrintLevel      Normal\n"
    "Print[ P_OrbEn ] 2   # 2 = print all orbital energies",
    "pal": "nprocs 8   # number of parallel cores",
    "paras": "R=1.05,1.13,8   # scan variable: start, end, number of steps\n"
    "# reference in the coordinate block as e.g. N 0 0 {R}",
    "plots": "dim1 45\ndim2 45\ndim3 45\nFormat Cube\n"
    'MO("orbital.cube",4,0);   # plot MO index 4',
    "rel": "method   ZORA   # DKH, ZORA, IORA, IORAmm, X2C\nModelPot 1,1,1,1",
    "rocis": "NRoots 6\nEWin = -5,5,-5,5   # donor/acceptor orbital energy window",
    "rr": "states 1,2,3,4,5\n"
    'HessName "ground_state.hess"\n'
    "ASAInput True\n"
    "# legacy module - the newer %esd block is recommended for new work",
    "scf": "MaxIter 500   # raise SCF iteration cap\n"
    "Guess   HCore   # HCore, Hueckel, PAtom, PModel, MORead",
    "symmetry": "SymThresh 0.01   # loosen point-group detection tolerance",
}


def format_block_display(block_name, content):
    """
    Build the %blockname/END-wrapped text shown for a block in
    List_add_job. No indentation is added (unlike the old display format)
    so the block name (first line) and END (last line) can be reliably
    stripped back off again when the user edits the item in place - see
    CalcSetupWindowORCA.on_block_item_edited.
    """
    return "%%%s\n%s\nEND" % (block_name, content)


class CalcSetupWindowORCA(QtWidgets.QMainWindow, Ui_SetupWindow):
    def __init__(self, parent, filepath):
        super(CalcSetupWindowORCA, self).__init__(parent)
        self.react = parent
        self.filepath = filepath
        self.settings = self.react.settings

        # Ensure we're using ORCA settings
        if self.settings.software != "ORCA":
            self.settings.software = "ORCA"

        # TODO this is temporary solution -> self.mol_obj should have a property state?
        self.state = self.react.get_current_state

        self.pymol = False
        self.pymol_animation = False
        if self.react.pymol:
            self.pymol = self.react.pymol
            self.pymol.pymol_cmd("set mouse_selection_mode, 0")
            self.pymol.pymol_cmd("hide cartoon")

            # Monitor clicks in pymol stdout_handler
            self.pymol.monitor_clicks()
            # Connect signals from pymol:
            self.pymol.atomsSelectedSignal.connect(self.pymol_atom_clicked)

            # Connect signal, update charge variable with the signal emitted
            self.pymol_charge = None
            self.pymol.overallChargeSignal.connect(self.handle_overall_charge)
            self.pymol.set_overall_charge()

        self.ui = Ui_SetupWindow()
        self.ui.setupUi(self)
        # TODO auto freeze function not working, hiding button for now
        self.ui.button_auto_freeze.hide()

        # TODO optimize this?:
        screen_size = self.screen().availableGeometry()
        window_size = self.geometry()
        self.move(int(((screen_size.width() - window_size.width()) / 2)) - 250, 100)

        self.setWindowTitle("REACT - Calculation setup")

        self.mol_obj = self.react.states[self.react.state_index].get_molecule_object(
            self.filepath
        )

        # Check if filetype is xyz (no charge information available)
        filetype = self.filepath.split(".")[-1].lower()
        self.is_xyz_file = filetype == "xyz"

        if self.is_xyz_file:
            self.react.append_text(
                "WARNING: XYZ files do not contain charge information. "
                "Please set the charge manually in the setup window."
            )
            self.charge = ""
        elif self.mol_obj.charge:
            self.charge = self.mol_obj.charge
        else:
            self.charge = None

        if self.mol_obj.multiplicity:
            self.multiplicity = self.mol_obj.multiplicity
        else:
            self.multiplicity = "1"
        self.filename = self.mol_obj.filename.split(".")[0]

        # we need to make a local copy of all info job-related stuff,
        # So that we dont change the attributes in Settings object
        self.functional = copy.deepcopy(self.settings.functional)
        self.basis = copy.deepcopy(self.settings.basis)
        self.job_type = self.settings.job_type
        self.opt_freq_combi = False
        self.job_options = copy.deepcopy(self.settings.job_options)
        self.blocks = copy.deepcopy(self.settings._orca_settings.get("blocks", {}))
        self.blocks_available = copy.deepcopy(
            self.settings._orca_settings.get("blocks_available", {})
        )
        self.solvents = ["CPCM", "SMD"]

        self.opt_freq_details = {"checked": False, "keywords": []}
        self.num_files = 1
        self.scan_bond_files = {}
        self.IRC_files = {}
        self.atom_bonds = {}
        self.Qbutton_scan_group = QtWidgets.QButtonGroup(self)
        self.Qbutton_scan_group.addButton(self.ui.radioButton_plus)
        self.Qbutton_scan_group.addButton(self.ui.radioButton_minus)
        self.Qbutton_scan_group.addButton(self.ui.radioButton_both)

        self.insert_model_atoms()
        self.fill_main_tab()

        self.ui.Button_add_job.clicked.connect(self.add_block)
        self.ui.List_add_job.itemChanged.connect(self.on_block_item_edited)
        self.ui.Button_del_job.clicked.connect(
            lambda: self.del_item_from_list(
                self.ui.List_add_job, self.job_options[self.job_type]
            )
        )

        self.ui.Button_add_job_2.clicked.connect(
            lambda: self.add_item_to_list(
                self.ui.LineEdit_add_job_2,
                self.ui.List_add_job_2,
                self.job_options[self.job_type],
            )
        )
        self.ui.Button_del_job_2.clicked.connect(
            lambda: self.del_item_from_list(self.ui.List_add_job_2, self.job_options)
        )
        self.ui.comboBox_basis1.currentIndexChanged.connect(self.update_basis1)
        self.ui.comboBox_funct.currentTextChanged.connect(self.update_functional)
        self.ui.comboBox_job_type.currentTextChanged.connect(self.update_job_details)
        self.ui.ComboBox_files.currentTextChanged.connect(self.update_preview_combobox)
        self.ui.button_close.clicked.connect(self.on_close)
        self.ui.button_write.clicked.connect(self.on_write)
        self.ui.button_add_freeze.clicked.connect(self.add_freeze_atoms)
        self.ui.button_delete_freeze.clicked.connect(self.remove_freeze_atoms)
        self.ui.button_auto_freeze.clicked.connect(self.auto_freeze_atoms)
        self.ui.button_add_scan.clicked.connect(self.add_scan_atoms)
        self.ui.button_delete_scan.clicked.connect(self.remove_scan_atoms)
        self.ui.lineEdit_filename.textChanged.connect(self.filename_update)
        self.ui.lineEdit_charge.textChanged.connect(self.update_charge)
        self.ui.lineEdit_multiplicity.textChanged.connect(self.update_multiplicity)
        self.ui.tabWidget.currentChanged.connect(self.update_preview)
        self.Qbutton_scan_group.buttonClicked.connect(self.on_scan_mode_changed)
        self.ui.spinbox_scan_pm.valueChanged.connect(
            lambda: self.on_spinbox_changed(self.ui.spinbox_scan_pm)
        )
        self.ui.spinbox_scan_increment.valueChanged.connect(
            lambda: self.on_spinbox_changed(self.ui.spinbox_scan_increment)
        )
        self.ui.button_invert.clicked.connect(self.on_invert_atoms)
        self.ui.checkBox_moveboth.clicked.connect(self.on_move_both_changed)

        self.ui.pushButton_create_mv.clicked.connect(self.add_move_atoms)
        self.ui.checkBox_moveboth_mv.clicked.connect(self.on_move_both_changed_mv)
        self.ui.spinbox_mv_bonds.valueChanged.connect(
            lambda: self.on_spinbox_changed(self.ui.spinbox_mv_bonds)
        )
        self.move_bond = False
        self.atoms_mv = list()
        self.ui.Button_save_new_geometry.clicked.connect(self.save_geometry_to_file)

        self.ui.list_model.itemSelectionChanged.connect(
            lambda: self.model_atom_clicked(self.ui.list_model)
        )
        self.ui.list_model.setSelectionMode(1)

        self.ui.list_model_mv.itemSelectionChanged.connect(
            lambda: self.model_atom_clicked(self.ui.list_model_mv)
        )
        self.ui.list_model_mv.setSelectionMode(2)

        self.ui.list_freeze_atoms.itemSelectionChanged.connect(self.freeze_list_clicked)

        self.atoms_to_select = 2
        self.ui.comboBox_freezetype.currentTextChanged.connect(
            self.change_selection_mode
        )

        self.ui.checkbox_eps.stateChanged.connect(self.update_eps)
        self.ui.comboBox_SCRF.currentTextChanged.connect(self.update_solvation_method)
        self.ui.comboBox.currentTextChanged.connect(self.update_block_preview)

        # Keep track of atoms selected (for multiple selection options)
        self.selected_indexes = list()

        # Pymol returns ordered atom id list, so we need to keep track of click order
        self.selected_ids = list()

        # key to get scan bond from list of bonds
        self.scan_bond_key = False

    @property
    def scan_bond(self):
        """
        Return AtomBond obj representing the bond selected from the scan list
        """
        try:
            return self.atom_bonds[self.scan_bond_key]
        except KeyError:
            pass

    def on_invert_atoms(self):
        if self.ui.list_model.count() < 1:
            return

        self.scan_bond.invert_atoms()
        self.ui.lineEdit_freeze.setText(str(self.scan_bond.atom1_idx))
        self.ui.lineEdit_move.setText(str(self.scan_bond.atom2_idx))
        self.update_scan()

    def on_move_both_changed(self):
        if self.ui.list_model.count() < 1:
            return

        if self.ui.checkBox_moveboth.isChecked():
            self.scan_bond.move_both = True
            disable = True
        else:
            self.scan_bond.move_both = False
            disable = False

        self.ui.button_invert.setDisabled(disable)
        self.ui.lineEdit_freeze.setDisabled(disable)
        self.ui.lineEdit_move.setDisabled(disable)
        self.update_scan()

    def on_move_both_changed_mv(self):
        if self.move_bond == False:
            return

        if self.ui.checkBox_moveboth_mv.isChecked():
            self.move_bond.move_both = True
        else:
            self.move_bond.move_both = False

        self.add_move_atoms()

    def update_eps(self):
        enabled = self.ui.checkbox_eps.isChecked()
        self.ui.lineEdit_eps.setEnabled(enabled)
        self.ui.lineEdit_refrac.setEnabled(enabled)

    def update_solvation_method(self):
        """
        Custom epsilon/refractive index and the cds_cpcm cavity/dispersion
        term are CPCM-specific (see ORCA 6 manual) - SMD uses its own,
        separate solvent-parameter library and reports its own CDS term.
        Disable (and uncheck, rather than leaving a stale checked-but-
        disabled state that would still affect the generated input) these
        controls whenever SMD is selected.
        """
        is_cpcm = self.ui.comboBox_SCRF.currentText() == "CPCM"

        if not is_cpcm:
            self.ui.checkbox_eps.setChecked(False)
            self.ui.checkbox_cds.setChecked(False)

        self.ui.checkbox_eps.setEnabled(is_cpcm)
        self.ui.checkbox_cds.setEnabled(is_cpcm)
        # lineEdit_eps/lineEdit_refrac stay governed by checkbox_eps (via
        # update_eps) - checkbox_eps being unchecked above already disables
        # them when switching away from CPCM.

    def update_block_preview(self):
        """
        Show a short, read-only example of the selected %block's most
        common/mandatory keywords, so the user has something concrete to
        start from rather than needing to already know ORCA's block
        syntax. Purely informational on its own - clicking Add (see
        add_block) is what actually copies this into a new, directly
        editable entry in List_add_job.
        """
        block_name = self.ui.comboBox.currentText()
        preview = BLOCK_KEYWORD_PREVIEWS.get(
            block_name,
            "No keyword preview available for this block - see the ORCA "
            "manual, %%%s section." % block_name,
        )
        self.ui.textEdit_block_preview.setPlainText(preview)

    def populate_block_list(self):
        """
        (Re)populate List_add_job from self.blocks with directly-editable
        items (see add_block/on_block_item_edited) - shared by
        fill_main_tab (initial load) and update_job_details (which clears
        and rebuilds the list on every job-type change).
        """
        self.ui.List_add_job.clear()
        for block_name, block_content in self.blocks.items():
            item = QtWidgets.QListWidgetItem(
                format_block_display(block_name, block_content)
            )
            item.setFlags(item.flags() | Qt.ItemIsEditable)
            item.setData(Qt.UserRole, block_name)
            self.ui.List_add_job.addItem(item)

    def add_block(self):
        """
        Add the currently-selected block type to List_add_job, seeded with
        its common-keywords preview content (blank if no verified preview
        exists - see BLOCK_KEYWORD_PREVIEWS) rather than the generic
        "no preview available" message, which isn't real ORCA syntax and
        shouldn't be inserted as if it were. The resulting list item is
        directly editable in place (see on_block_item_edited) so the user
        adjusts it there instead of via a separate input field.
        If this block type is already in the list, update that entry in
        place instead of adding a duplicate.
        """
        block_name = self.ui.comboBox.currentText()
        if not block_name:
            return

        content = BLOCK_KEYWORD_PREVIEWS.get(block_name, "")
        self.blocks[block_name] = content
        display_text = format_block_display(block_name, content)

        for i in range(self.ui.List_add_job.count()):
            item = self.ui.List_add_job.item(i)
            if item.data(Qt.UserRole) == block_name:
                item.setText(display_text)
                return

        item = QtWidgets.QListWidgetItem(display_text)
        item.setFlags(item.flags() | Qt.ItemIsEditable)
        item.setData(Qt.UserRole, block_name)
        self.ui.List_add_job.addItem(item)

    def on_block_item_edited(self, item):
        """
        Keep self.blocks in sync when the user edits a block's entry in
        List_add_job directly, including renaming it by editing its
        "%blockname" first line. The block name is re-parsed from the
        edited text (rather than trusting the item's stored UserRole,
        which only tracks what it was before this edit) so a rename is
        picked up correctly; the stored UserRole is then refreshed to the
        new name for the next edit.
        """
        lines = item.text().split("\n")
        while lines and lines[-1] == "":
            lines.pop()

        if not lines or not lines[0].startswith("%"):
            # Malformed (user deleted the %blockname line) - leave
            # self.blocks alone rather than guess what they meant.
            return

        new_block_name = lines[0][1:].strip()
        if not new_block_name:
            return

        if lines[-1].strip().upper() == "END":
            content = "\n".join(lines[1:-1])
        else:
            # User deleted the END line - keep the rest as content anyway,
            # format_block_display will add a fresh END back when this
            # block is next displayed/regenerated.
            content = "\n".join(lines[1:])

        old_block_name = item.data(Qt.UserRole)
        if old_block_name and old_block_name != new_block_name:
            self.blocks.pop(old_block_name, None)

        self.blocks[new_block_name] = content
        item.setData(Qt.UserRole, new_block_name)

    def change_selection_mode(self):
        if self.ui.comboBox_freezetype.currentText() == "Atoms":
            self.ui.list_model.setSelectionMode(1)
        else:
            self.ui.list_model.setSelectionMode(2)

        select = {"Atom": 1, "Bond": 2, "Angle": 3, "Dihedral": 4}
        self.atoms_to_select = select[self.ui.comboBox_freezetype.currentText()]

    def enable_scan(self, enable=True):
        self.ui.button_add_scan.setEnabled(enable)
        self.ui.button_delete_scan.setEnabled(enable)
        self.ui.spinbox_radius.setEnabled(enable)
        self.ui.spinbox_scan_pm.setEnabled(enable)
        self.ui.spinbox_scan_increment.setEnabled(enable)

    @pyqtSlot(list)
    def pymol_atom_clicked(self, ids):
        """
        Get atoms selected in pymol and select deselect in model atom listwidget
        """
        if not ids or None in ids:
            return

        atom_list = None
        if self.ui.tabWidget.currentIndex() == 1:
            atom_list = self.ui.list_model_mv
        elif self.ui.tabWidget.currentIndex() == 2:
            atom_list = self.ui.list_model
        
        # If not on a relevant tab, ignore the signal
        if atom_list is None:
            return

        ids = [int(x) - 1 for x in ids]

        # Effective when len(ids) > len(self.selected_ids)
        unsele = list(set(self.selected_ids) - set(ids))
        try:
            for id in unsele:
                # Check if id is within range before accessing
                if id < 0 or id >= atom_list.count():
                    continue
                item = atom_list.item(id)
                if item is not None:
                    item.setSelected(False)
                    self.selected_ids.pop(self.selected_ids.index(id))
        except (UnboundLocalError, ValueError):
            pass

        # Effective when len(ids) < len(self.selected_ids)
        new_select = list(set(ids) - set(self.selected_ids))
        for id in new_select:
            self.selected_ids.append(id)
            if len(self.selected_ids) > self.atoms_to_select:
                self.selected_ids.pop(0)
            # Check if id is within range before accessing
            if id >= 0 and id < atom_list.count():
                item = atom_list.item(id)
                if item is not None:
                    item.setSelected(True)

    def model_atom_clicked(self, atom_list):
        """
        When selection i atom list is changed, update, and communicate with pymol
        """
        if self.pymol_animation:
            self.stop_pymol_animation()

        sele = atom_list.selectedIndexes()
        while len(sele) > self.atoms_to_select:
            if len(self.selected_indexes) > 0:
                atom_list.item(self.selected_indexes.pop(0).row()).setSelected(False)
            else:
                atom_list.item(sele.pop(0).row()).setSelected(False)
            sele = atom_list.selectedIndexes()

        self.selected_indexes = sele

        atoms = list()
        coordinates = list()
        for i in self.selected_indexes:
            atoms.append(self.mol_obj.molecule[i.row() + 1].atom_index)
            coordinates.append(self.mol_obj.molecule[i.row() + 1].coordinate)

        if len(atoms) == 2:
            # self.enable_scan(enable=True)
            r = atom_distance(coordinates[0], coordinates[1])
            self.ui.spinbox_radius.setValue(r)
            self.ui.spinbox_curr_mv.setValue(r)
        # else:
        #    self.enable_scan(enable=False)

        if self.pymol:
            self.update_pymol_selection(atoms=atoms)

    def parse_atom_indices(self, text):
        indexes = []
        if ":" in text:
            try:
                # Format: "Type: 0, 1, 2"
                tokens = text.split(":", 1)[1].split(",")
                for x in tokens:
                    if x.strip():
                        indexes.append(int(x.strip()))
            except ValueError:
                pass
        elif "(ORCA" in text:
            tokens = text.split()
            for i, token in enumerate(tokens):
                if token == "(ORCA" and i > 0:
                    try:
                        indexes.append(int(tokens[i - 1]))
                    except ValueError:
                        pass
        else:
            tokens = text.split()
            if len(tokens) > 2:
                for x in tokens[1:-1]:
                    try:
                        indexes.append(int(x))
                    except ValueError:
                        pass
        return indexes

    def freeze_list_clicked(self):
        """
        When entry in "Atoms to freeze" is clicked, update to selected in "Atoms in model" list and pymol
        """
        if self.pymol_animation:
            self.stop_pymol_animation()

        try:
            text = self.ui.list_freeze_atoms.currentItem().text()
            # Stored indices are 0-based
            indexes = self.parse_atom_indices(text)
        except AttributeError:
            return
        freeze_type = {1: "Atom", 2: "Bond", 3: "Angle", 4: "Dihedral"}
        if len(indexes) in freeze_type:
            self.ui.comboBox_freezetype.setCurrentText(freeze_type[len(indexes)])

        self.selected_indexes = indexes
        self.ui.list_model.clearSelection()

        for index in indexes:
            self.ui.list_model.item(index).setSelected(True)

    def update_pymol_selection(self, atoms):
        filetype = self.mol_obj.filepath.split(".")[-1]
        object_name = self.mol_obj.molecule_name + "_" + filetype
        group = "state_%d" % self.react.get_current_state
        self.pymol.set_selection(
            atoms=atoms, sele_name="sele", object_name=object_name, group=group
        )

    def add_freeze_atoms(self):
        """
        Adds selected atoms to freeze section - append F to end of str for freeze
        """
        if len(self.selected_indexes) < 1:
            return
        type = "C"
        g_cmd = {"Atom": "C", "Bond": "B", "Angle": "A", "Dihedral": "D"}
        atoms = []
        for i in self.selected_indexes:
            atomnr = self.mol_obj.molecule[i.row() + 1].atom_index
            if self.pymol:
                self.pymol_spheres(atomnr)
            atoms.append(f"{atomnr - 1}")

        atoms_str = ", ".join(atoms)
        constraint_str = f"{self.ui.comboBox_freezetype.currentText()}: {atoms_str}"

        for i in range(self.ui.list_freeze_atoms.count()):
            if self.ui.list_freeze_atoms.item(i).text() == constraint_str:
                return

        self.ui.list_freeze_atoms.insertItem(0, constraint_str)

    def remove_freeze_atoms(self):
        """
        Removes selected freeze atoms/bonds/angles/torsions
        """
        # Get rows for selected to delete:
        to_del = [
            self.ui.list_freeze_atoms.row(x)
            for x in self.ui.list_freeze_atoms.selectedItems()
        ]
        for row in to_del:
            if self.pymol:
                text = self.ui.list_freeze_atoms.item(row).text()
                for i in self.parse_atom_indices(text):
                    # convert back to 1-based for pymol
                    self.pymol_spheres(atom_nr=i + 1, hide=True)

            self.ui.list_freeze_atoms.takeItem(row)

    def on_spinbox_changed(self, spinbox):
        try:
            if spinbox.value() == 0:
                return
            if spinbox == self.ui.spinbox_scan_pm:
                self.scan_bond.scan_dist = spinbox.value()
            elif spinbox == self.ui.spinbox_scan_increment:
                self.scan_bond.step_size = spinbox.value()
            elif spinbox == self.ui.spinbox_mv_bonds:
                self.move_bond.scan_dist = spinbox.value()
        except AttributeError:
            return

        self.update_scan()

    def add_scan_atoms(self):
        """
        First, gather information from the UI. Then, pass all information
        to a new instance of AtomBond object.
        """

        if len(self.selected_indexes) != 2:
            self.react.append_text("Scan only works when only 2 atoms are selected...")
            return

        if self.scan_bond_key:
            try:
                del self.atom_bonds[self.scan_bond_key]
            except:
                pass

        freeze = ""
        atoms = list()
        bond_size = self.ui.spinbox_radius.value()
        scan_size = self.ui.spinbox_scan_pm.value()
        scan_increment = self.ui.spinbox_scan_increment.value()

        for i in self.selected_indexes:
            atoms.append(self.mol_obj.molecule[i.row() + 1].atom_index)
            freeze += f"{self.mol_obj.molecule[i.row() + 1].atom_name} ({self.mol_obj.molecule[i.row() + 1].atom_index})    "

        if self.ui.radioButton_plus.isChecked():
            scan_mode = "+"
        elif self.ui.radioButton_minus.isChecked():
            scan_mode = "-"
        else:
            scan_mode = "+/-"

        self.atom_bonds[freeze] = AtomBond(
            self.mol_obj.formatted_xyz,
            atoms[0],
            atoms[1],
            bond_size,
            scan_size,
            scan_increment,
            self.ui.checkBox_moveboth.isChecked(),
            scan_mode,
        )

        self.scan_bond_key = freeze

        # list_item = QtWidgets.QListWidgetItem(freeze)
        # self.ui.list_scan_bonds.addItem(list_item)
        # self.ui.list_scan_bonds.setCurrentItem(list_item)
        self.ui.lineEdit_freeze.setText(str(atoms[0]))
        self.ui.lineEdit_move.setText(str(atoms[1]))

        self.update_scan()

    def add_move_atoms(self):
        """
        First, gather information from the UI. Then, pass all information
        to a new instance of AtomBond object.
        """

        bond_size = self.ui.spinbox_curr_mv.value()
        scan_size = self.ui.spinbox_mv_bonds.value()
        scan_increment = self.ui.spinbox_scan_increment.value()

        self.atoms_mv.clear()

        for i in self.selected_indexes:
            self.atoms_mv.append(self.mol_obj.molecule[i.row() + 1].atom_index)

        if self.ui.radioButton_plus_mv.isChecked():
            scan_mode = "+"
        else:
            scan_mode = "-"

        try:
            self.move_bond = AtomBond(
                self.mol_obj.formatted_xyz,
                self.atoms_mv[0],
                self.atoms_mv[1],
                bond_size,
                scan_size,
                scan_increment,
                self.ui.checkBox_moveboth_mv.isChecked(),
                scan_mode,
            )
        except:
            pass

        all_xyz = self.move_bond.scan_new_coordinates
        last_xyz = all_xyz.popitem()

        content = str(len(last_xyz[1])) + "\n" + last_xyz[0] + "\n"
        for line in last_xyz[1]:
            content += line + "\n"

        filepath = self.settings.workdir + "/.move_tmp.xyz"
        try:
            remove(filepath)
        except:
            pass

        with open(filepath, "w+") as f:
            f.write(content)

        self.react.pymol.load_structure(filepath, delete_after=True)
        self.react.pymol.pymol_cmd(
            "group state_%d, %s" % (1, filepath.split("/")[-1].split(".")[0])
        )

        self.react.pymol.set_default_rep()
        self.react.pymol.pymol_cmd(
            "enable state_%d and %s" % (1, filepath.split("/")[-1].split(".")[0])
        )
        # self.react.file_to_pymol(self.settings.workdir + '/.move_tmp.xyz', state=1, set_defaults=True)

        # if self.pymol:
        #    self.react.add_file(self.settings.workdir + '/.move_tmp.xyz')

    def save_geometry_to_file(self):
        bond_size = self.ui.spinbox_curr_mv.value()
        scan_size = self.ui.spinbox_mv_bonds.value()
        scan_increment = self.ui.spinbox_scan_increment.value()

        for i in self.selected_indexes:
            self.atoms_mv.append(self.mol_obj.molecule[i.row() + 1].atom_index)

        if self.ui.radioButton_plus_mv.isChecked():
            scan_mode = "+"
        else:
            scan_mode = "-"

        try:
            self.move_bond = AtomBond(
                self.mol_obj.formatted_xyz,
                self.atoms_mv[0],
                self.atoms_mv[1],
                bond_size,
                scan_size,
                scan_increment,
                self.ui.checkBox_moveboth_mv.isChecked(),
                scan_mode,
            )
        except IndentationError:
            pass

        all_xyz = self.move_bond.scan_new_coordinates
        last_xyz = all_xyz.popitem()

        content = str(len(last_xyz[1])) + "\n" + last_xyz[0] + "\n"
        for line in last_xyz[1]:
            content += line + "\n"

        filepath = self.settings.workdir + f"/{last_xyz[0]}.xyz"

        with open(filepath, "w+") as f:
            f.write(content)

        self.react.add_file(filepath)

    def update_scan(self):
        """
        Remove all old files before making new ones
        """
        if not path.isdir(f"{self.settings.workdir}/.scan_temp"):
            mkdir(f"{self.settings.workdir}/.scan_temp")
        else:
            tempfiles = glob.glob(f"{self.settings.workdir}/.scan_temp/*")
            for f in tempfiles:
                remove(f)

        try:
            self.scan_bond.write_xyzfiles(self.settings.workdir + "/.scan_temp/")
            if self.pymol:
                # read files to pymol obj
                self.anmiate_bond_pymol()
        except AttributeError:
            pass

    def update_move(self):
        """
        Remove all old files before making new ones
        """
        # if path.isfile(f"{self.settings.workdir}/.move_temp"):
        #    self.react.remove_file(f"{self.settings.workdir}/.move_temp")
        #    remove(f"{self.settings.workdir}/.move_temp")
        #    self.ui.list_model_mv.clearSelection()

        self.add_move_atoms()

    def anmiate_bond_pymol(self):
        # It is not possible to delete states, so we actually need to delete the mol_obj and load it again :/
        # name_split = self.mol_obj.filepath.split("/")[-1].split(".")
        # name = name_split[0] + "_" + name_split[-1]
        # self.pymol.pymol_cmd(f"delete {name} and state_{self.state}")
        # self.react.file_to_pymol(filepath=self.mol_obj.filepath, state=1, set_defaults=True)
        state = 1
        self.pymol.pymol_cmd(f"delete scan")
        for f in sorted(os.listdir(f"{self.settings.workdir}/.scan_temp")):
            state += 1
            _mol = f"{self.settings.workdir}/.scan_temp/{f}"
            self.pymol.pymol_cmd(f"load {_mol},  scan, {state}")
        self.pymol.set_default_rep()
        self.pymol.pymol_cmd("set movie_fps, 10")
        self.pymol_animation = True
        self.pymol.pymol_cmd("mplay")

    def stop_pymol_animation(self):
        self.pymol.pymol_cmd("mstop")
        self.pymol.pymol_cmd(f"delete scan")
        # self.react.file_to_pymol(filepath=self.mol_obj.filepath, state=1, set_defaults=True)
        self.pymol_animation = False

    def remove_scan_atoms(self):
        try:
            del self.atom_bonds[self.scan_bond_key]
            self.ui.lineEdit_move.clear()
            self.ui.lineEdit_freeze.clear()
            self.ui.spinbox_radius.clear()
        except:
            pass
        # self.update_scan()

    # def scan_list_clicked(self):
    #     if not hasattr(self, "scan_bond"):
    #         return
    #     self.ui.list_model_scan.clearSelection()

    #     self.ui.lineEdit_freeze.setText(str(self.scan_bond.atom1_idx))
    #     self.ui.lineEdit_move.setText(str(self.scan_bond.atom2_idx))

    #     self.ui.checkBox_moveboth.setChecked(self.scan_bond.move_both)
    #     self.ui.spinbox_radius.setValue(self.scan_bond.bond_dist)
    #     self.ui.spinbox_scan_pm.setValue(self.scan_bond.scan_dist)
    #     self.ui.spinbox_scan_increment.setValue(self.scan_bond.step_size)
    #     if self.scan_bond.scan_mode == '+':
    #         self.ui.radioButton_plus.setChecked(True)
    #     elif self.scan_bond.scan_mode == '-':
    #         self.ui.radioButton_plus.setChecked(True)
    #     else:
    #         self.ui.radioButton_both.setChecked(True)

    #     for i in [self.scan_bond.atom1_idx - 1, self.scan_bond.atom2_idx - 1]:
    #         self.ui.list_model_scan.item(i).setSelected(True)

    #     #self.scan_bond.write_xyzfiles()

    def pymol_spheres(self, atom_nr, hide=False):
        """
        Indicate withe spheres atoms to be frozen / modredundant
        """
        group = "state_%d" % self.react.get_current_state
        _cmd = "show"
        mol_obj_name = (
            self.mol_obj.filename.split(".")[0]
            + "_"
            + self.mol_obj.filepath.split(".")[-1]
        )
        if hide:
            _cmd = "hide"
        self.pymol.pymol_cmd(
            f"{_cmd} spheres, id {atom_nr} and {group} and {mol_obj_name}"
        )
        if not hide:
            self.pymol.pymol_cmd(
                f"set sphere_scale, 0.3, id {atom_nr} and {group} and {mol_obj_name}"
            )

    def update_job_details(self):
        """
        Activated when job option combobox is updated. Will fill Qlist,
        Enable multiple files in case IRC is chosen, and hide Opt+freq
        buttons whenever job != opt.
        :return:
        """
        self.IRC_files.clear()
        self.ui.ComboBox_files.blockSignals(True)
        self.ui.ComboBox_files.clear()
        self.ui.ComboBox_files.blockSignals(False)

        self.job_type = self.ui.comboBox_job_type.currentText()

        if self.job_type in ["Freq"]:
            self.ui.checkbox_freq.setHidden(True)
        else:
            self.ui.checkbox_freq.setHidden(False)

        self.ui.List_add_job_2.clear()

        self.ui.List_add_job_2.addItems(self.job_options[self.job_type])
        self.populate_block_list()

        # Automatically select 'neb' block for NEB-TS job type
        if self.job_type == "NEB-TS":
            neb_index = self.ui.comboBox.findText("neb")
            if neb_index >= 0:
                self.ui.comboBox.setCurrentIndex(neb_index)

    def is_composite_method(self):
        """Check if current functional is a composite method"""
        functional_options = self.settings.functional_options
        # Ensure functional_options is a dict (not list from old settings)
        if not isinstance(functional_options, dict):
            return False
        if self.functional in functional_options:
            func_info = functional_options[self.functional]
            if isinstance(func_info, dict):
                return func_info.get("Composite method", False)
        return False

    def update_functional(self):
        self.functional = self.ui.comboBox_funct.currentText()
        # Disable basis selection for composite methods
        if self.is_composite_method():
            self.ui.comboBox_basis1.setEnabled(False)
            # Make it visually grayed out
            self.ui.comboBox_basis1.setStyleSheet(
                "QComboBox { color: gray; background-color: rgb(40, 40, 40); }"
            )
        else:
            self.ui.comboBox_basis1.setEnabled(True)
            # Restore normal appearance
            self.ui.comboBox_basis1.setStyleSheet("")

    def fill_main_tab(self):
        """
        Fill all widgets in main tab.
        ORCA doesn't use basis set modifiers - they're built into the basis name.
        :return:
        """
        self.ui.lineEdit_filename.setText(self.filename)

        # Get software-specific options
        functional_options = self.settings.functional_options
        basis_options = self.settings.basis_options

        self.ui.comboBox_funct.addItems(list(functional_options.keys()))
        self.ui.comboBox_basis1.addItems(basis_options)

        self.ui.comboBox_job_type.addItems(self.job_options)
        self.populate_block_list()
        self.ui.ComboBox_files.addItem(self.filename)

        self.ui.comboBox.addItems([x for x in self.blocks_available.keys()])
        self.update_block_preview()

        self.ui.comboBox_job_type.setCurrentText(self.settings.job_type)
        self.ui.comboBox_funct.setCurrentText(self.functional)
        self.ui.comboBox_basis1.setCurrentText(self.basis)

        self.ui.lineEdit_charge.setText(self.charge)
        self.ui.lineEdit_multiplicity.setText(self.multiplicity)

        self.ui.comboBox_SCRF.addItems(self.solvents)
        self.ui.checkbox_eps.setChecked(False)
        self.ui.lineEdit_eps.setEnabled(False)
        self.ui.lineEdit_refrac.setEnabled(False)
        self.ui.checkbox_draco.setChecked(False)
        self.ui.checkbox_cds.setChecked(False)
        self.update_solvation_method()

        self.update_job_details()

        # Ensure basis combobox state is set correctly based on functional
        self.update_functional()

    def on_scan_mode_changed(self):
        if not hasattr(self, "scan_bond"):
            return
        scan_mode_map = {-2: "+", -3: "-", -4: "+/-"}
        key = self.Qbutton_scan_group.checkedId()
        self.scan_bond.scan_mode = scan_mode_map[key]
        self.update_scan()

    def update_preview(self):
        """
        Update preview tab. Combobox is cleared and loaded again everytime.
        In case of multiple files, add correct num of items to combobox,
        but make only file content for the last file in combobox, and load
        this to text preview box.
        """

        if self.ui.tabWidget.currentIndex() == 1:
            self.atoms_to_select = 2

        # check if preview tab is selected. if not, return
        if not self.ui.tabWidget.currentIndex() == 3:
            return

        self.ui.ComboBox_files.blockSignals(True)
        self.ui.ComboBox_files.clear()
        self.ui.ComboBox_files.blockSignals(False)

        if self.atom_bonds:
            for bond, bond_obj in self.atom_bonds.items():
                for filename, xyz in bond_obj.scan_new_coordinates.items():
                    self.ui.ComboBox_files.addItem(filename)

                self.ui.ComboBox_files.setCurrentText(filename)
                file_content = self.make_input_content(
                    filename=filename, xyz=xyz, bond_obj=bond_obj
                )
        elif self.IRC_files:
            for filename, keywords in self.IRC_files.items():
                self.ui.ComboBox_files.addItem(filename)

            self.ui.ComboBox_files.setCurrentText(filename)
            file_content = self.make_input_content(
                filename=filename, extra_job_keywords=keywords
            )
        else:
            file_content = self.make_input_content(filename=self.filename)
            self.ui.ComboBox_files.addItem(self.filename)

        self.ui.text_preview.setPlainText(file_content)

    def update_preview_combobox(self):
        """
        This function is called everytime any change to combobox is made.
        Program crashes if runs with only one file in comobox, thus the
        if statement is neccesarry!
        """
        if self.ui.ComboBox_files.count() > 1:
            filename = self.ui.ComboBox_files.currentText()
            try:
                file_content = self.make_input_content(
                    filename=filename, extra_job_keywords=self.IRC_files[filename]
                )
                self.ui.text_preview.setPlainText(file_content)
                return
            except KeyError:
                pass

            for bond, bond_obj in self.atom_bonds.items():
                for filename_scan, xyz in bond_obj.scan_new_coordinates.items():
                    if filename_scan == filename:
                        file_content = self.make_input_content(
                            filename=filename, xyz=xyz, bond_obj=bond_obj
                        )
                        self.ui.text_preview.setPlainText(file_content)

    def update_charge(self):
        self.charge = self.ui.lineEdit_charge.text()

    def update_multiplicity(self):
        self.multiplicity = self.ui.lineEdit_multiplicity.text()

    def make_files(self):
        """
        Writes one or multiple inputfiles (found in self.multiplefiles)
        New files are loaded into react again using the exsisting add_files()
        function in REACT.py
        """
        files = []
        if self.IRC_files:
            for filename, keywords in self.IRC_files.items():
                content = self.make_input_content(
                    filename=filename, extra_job_keywords=keywords
                )
                filepath = self._make_file(filename, content)
                files.append(filepath)
        elif self.atom_bonds:
            for bond, bond_obj in self.atom_bonds.items():
                for filename_scan, xyz in bond_obj.scan_new_coordinates.items():
                    content = self.make_input_content(
                        filename=filename_scan, xyz=xyz, bond_obj=bond_obj
                    )
                    filepath = self._make_file(filename_scan, content)
                    files.append(filepath)
        else:
            content = self.make_input_content(filename=self.filename)
            filepath = self._make_file(self.filename, content)
            files.append(filepath)

        # Add files separate to avoid issues relating to multithreading prosess in REACT.py
        for filepath in files:
            if filepath == None:
                # Something happend (maybe user changed their mind and closed QfileDialog)
                return
            self.react.append_text(f"{filepath.split('/')[-1]} saved to: {filepath}")
            if self.ui.checkBox_cp_to_reactmain.isChecked():
                self.react.add_file(filepath)

    def make_orca_input_content(
        self, filename, extra_job_keywords=False, xyz=False, bond_obj=False
    ):
        """
        Make content (not file) for one ORCA inputfile.
        :return: str
        ORCA simple input format:
        ! functional basis job_keywords
        %blocks
        * xyz charge multiplicity
        coordinates
        *
        """
        # Simple input line
        job_keywords = []
        blocks_str = ""

        job_type = self.job_type
        if self.ui.checkbox_freq.isChecked():
            job_keywords.append("Freq")

        if extra_job_keywords:
            job_keywords.extend(extra_job_keywords)

        # Add job type specific keywords
        if self.job_type in self.job_options:
            job_keywords.extend(self.job_options[self.job_type])

        if self.ui.checkbox_SCRF.isChecked():
            scrf_method = self.ui.comboBox_SCRF.currentText()
            solvent = self.ui.lineEdit_solvent.text()
            # Custom epsilon/refrac and cds_cpcm are CPCM-specific (ORCA 6
            # manual) - SMD uses its own solvent-parameter library and
            # reports its own, separate CDS term. The UI already disables
            # these controls for SMD; this guard is just defensive.
            is_cpcm = scrf_method == "CPCM"

            cpcm_block_lines = []

            if is_cpcm and self.ui.checkbox_eps.isChecked():
                cpcm_block_lines.append(f"  epsilon {self.ui.lineEdit_eps.text()}")
                cpcm_block_lines.append(f"  refrac {self.ui.lineEdit_refrac.text()}")
                solvent = None

            if is_cpcm and self.ui.checkbox_cds.isChecked():
                # Since ORCA 6, plain CPCM no longer adds the cavity/
                # dispersion (cav+disp) free-energy term by default -
                # cds_cpcm 2 restores it (ORCA 6 manual, solvationmodels).
                cpcm_block_lines.append("  cds_cpcm 2")

            if cpcm_block_lines:
                blocks_str += "\n%cpcm\n" + "\n".join(cpcm_block_lines) + "\nEND\n"

            if solvent:
                job_keywords.append(f"{scrf_method}({solvent})")
            else:
                job_keywords.append(f"{scrf_method}")

            if self.ui.checkbox_draco.isChecked():
                job_keywords.append("DRACO")

        # Build simple input line
        # Composite methods don't need basis set specification
        if self.is_composite_method():
            simple_input = f"! {self.functional} {job_type} {' '.join(job_keywords)}"
        else:
            simple_input = (
                f"! {self.functional} {self.basis} {job_type} {' '.join(job_keywords)}"
            )

        if self.blocks:
            blocks_str += "\n" + ("\n").join(
                [
                    f"%{block_name}\n  {block_content}\nEND\n"
                    for block_name, block_content in self.blocks.items()
                ]
            )

        # Molecule specification
        charge_value = self.charge
        if not self.charge or self.charge == "":
            charge_value = "INSERT_CHARGE_HERE"
            if self.is_xyz_file:
                self.react.append_text(
                    "WARNING: XYZ file - no charge information. Please replace INSERT_CHARGE_HERE in the input file."
                )
            else:
                self.react.append_text(
                    "WARNING: No charge is found the system. Please replace INSERT_CHARGE_HERE in the input file."
                )

        if xyz:
            coords_str = "\n".join(xyz)
        else:
            coords_str = "\n".join(self.mol_obj.formatted_xyz)

        molecule_str = f"* xyz {charge_value} {self.multiplicity}\n{coords_str}\n*"

        # Geometry constraints (modredundant equivalent)
        constraints_str = ""
        if self.ui.list_freeze_atoms.count() > 0:
            constraints_str = "\n%geom\n"
            constraints_str += "  Constraints\n"
            for i in range(self.ui.list_freeze_atoms.count()):
                item = self.ui.list_freeze_atoms.item(i).text()
                # Parse format: "Constraint_type: 1, 2, 3"
                if ":" in item:
                    constraint_type_map = {
                        "Atom": "C",
                        "Bond": "B",
                        "Angle": "A",
                        "Dihedral": "D",
                    }
                    constraint_name, atoms_part = item.split(":", 1)
                    constraint_type = constraint_type_map.get(
                        constraint_name.strip(), "C"
                    )

                    # Extract atom numbers before parentheses
                    atom_indices = []
                    for atom_str in atoms_part.split(","):
                        atom_str = atom_str.strip()

                        atom_indices.append(atom_str)

                    if atom_indices:
                        constraints_str += (
                            f"    {{ {constraint_type} {' '.join(atom_indices)} C }}\n"
                        )
            constraints_str += "  END\nEND\n"

        return f"{simple_input}\n{blocks_str}\n{molecule_str}\n{constraints_str}\n"

    def _make_file(self, filename, file_content):
        """
        Private function. Makes sure that no files are overwritten, and
        makes an unique filepath by adding _1, _2, _3 etc until an unique
        filepath is found.
        """
        i = 0
        # ORCA uses .inp extension
        ext = ".inp"
        new_filepath = self.react.settings.workdir + "/" + filename + ext

        if path.isfile(new_filepath) == True:
            new_filepath, filter_ = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Filename already exists, please select a different filename",
                self.react.settings.workdir,
                "input files (*.com *.inp)",
            )
        if new_filepath == "":
            return
        with open(new_filepath, "w+") as f:
            f.write(file_content)
            f.write("\n")

        return new_filepath

    def make_input_content(
        self, filename, extra_job_keywords=False, xyz=False, bond_obj=False
    ):
        """
        Make content (not file) for one inputfile (Gaussian or ORCA).
        Routes to appropriate method based on software selection.
        :return: str
        """
        return self.make_orca_input_content(filename, extra_job_keywords, xyz, bond_obj)

    def remove_extra_newline(self, input_string):
        """
        To remedy some strange bug in make_input_content..
        """

        # Split the string into lines
        lines = input_string.split("\n")

        # Remove any empty lines
        lines = [line for line in lines if line]

        # Join the lines back together with newline characters
        output_string = "\n".join(lines)

        return output_string

    def route_checkboxes_update(self, checkbox, lineEdit):
        if checkbox.isChecked():
            lineEdit.setEnabled(True)
        else:
            lineEdit.setEnabled(False)

    def update_basis1(self):
        """
        Updates self.basis according to basis selected by user.
        ORCA doesn't use basis set modifiers (diff, pol1, pol2).
        """
        self.basis = self.ui.comboBox_basis1.currentText()

    def filename_update(self):
        self.filename = self.ui.lineEdit_filename.text()
        # ORCA doesn't use .chk files like Gaussian

    def insert_model_atoms(self):
        """
        :return:
        """

        if "pdb" in self.filepath.split(".")[-1]:
            atoms = self.mol_obj.formatted_pdb
        else:
            self.ui.button_auto_freeze.setEnabled(False)
            atoms = self.mol_obj.formatted_xyz

        for atom_list in [self.ui.list_model, self.ui.list_model_mv]:
            for i in range(len(atoms)):
                atom_list.insertItem(i, f"{i}  {atoms[i]}")

    def auto_freeze_atoms(self):
        """
        Goes through PDB atoms and tries to figure out base on pdb atom names and distances what atoms
        are terminal/chopped region atoms that should be frozen.
        """
        expected = ["CA", "C", "H", "N", "O"]
        atoms = self.mol_obj.atoms
        for atom in atoms:
            pdb_name = atom.get_pdb_atom_name
            if pdb_name in ["C", "N"]:
                atoms.pop(0)
                for next_atom in atoms:
                    radius = atom_distance(
                        atom.get_coordinate, next_atom.get_coordinate
                    )
                    if radius < 1.7:
                        if next_atom.get_pdb_atom_name not in expected:
                            atom_nr = next_atom.get_atom_index
                            self.ui.list_freeze_atoms.insertItem(0, f"X {atom_nr} F")
                            if self.pymol:
                                self.pymol_spheres(atom_nr)

    def add_item_to_list(self, Qtextinput, Qlist, job_list):
        """
        :param Qtextinput: QLineEdit
        :param Qlist: QListWidget
        :param job_list: list to store item
        Adds the text input from user (past to Qtextinput) to correct
        QlistWidget and append item to list.
        (%Block entries are handled separately - see add_block/
        on_block_item_edited - since they're keyed by block name in a
        dict rather than appended to a plain list.)
        """
        user_input = Qtextinput.text()

        if not user_input:
            return

        if user_input not in job_list:
            job_list.append(user_input)
            Qlist.addItem(user_input)

        Qtextinput.clear()

    def del_item_from_list(self, Qlist, job_list):
        """
        :param Qlist: QListWidget
        :param job_list: list to delete from
        Removes item in QlistWdiget and updates self.settings accordingly.
        """
        try:
            item_text = Qlist.currentItem().text()

            # For blocks list (List_add_job), extract block name and delete from self.blocks
            if Qlist == self.ui.List_add_job:
                # Extract block name from formatted string like "%blockname\n    content\nEND\n"
                block_name = item_text.split("\n")[0][1:]  # Remove '%' prefix
                if block_name in self.blocks:
                    del self.blocks[block_name]
            # For job options list (List_add_job_2), delete from job_options for current job type
            elif Qlist == self.ui.List_add_job_2:
                if (
                    self.job_type in self.job_options
                    and item_text in self.job_options[self.job_type]
                ):
                    self.job_options[self.job_type].remove(item_text)
            # For other lists, delete from the passed job_list
            elif item_text in job_list:
                job_list.remove(item_text)

            Qlist.takeItem(Qlist.currentRow())
        except (AttributeError, IndexError):
            pass

    def del_tempfiles(self):
        try:
            shutil.rmtree(self.react.react_path + "/.scan_temp")
        except OSError as e:
            pass

    def on_close(self):
        self.del_tempfiles()
        self.close()

    def on_write(self):
        """
        Write all data to Inputfile object.
        # TODO how to return this object to react after closing the window?
        # TODO reply --> just open the .com file in the current state (tab) of REACT after write?
        # TODO on_write should probably just write whatever is in the preview window to a file.
        """
        # self.mol_obj.filename = self.ui.lineEdit_filename.text()
        # self.mol_obj.job_type = self.ui.comboBox_job_type.currentText()
        # self.mol_obj.basis = self.ui.comboBox_basis1.currentText()
        # self.mol_obj.basis_diff = self.ui.comboBox_basis2.currentText()
        # self.mol_obj.basis_pol1 = self.ui.comboBox_basis3.currentText()
        # self.mol_obj.basis_pol2 = self.ui.comboBox_basis4.currentText()
        # self.react.states[self.react.state_index].add_instance(self.mol_obj)

        self.make_files()
        self.del_tempfiles()
        # self.close()

    def closeEvent(self, event):
        if self.pymol:
            self.pymol.pymol_cmd("set mouse_selection_mode, 1")
            self.pymol.pymol_cmd("hide spheres,")
            self.pymol.unmonitor_clicks()
        if path.isdir(f"{self.settings.workdir}/.scan_temp"):
            shutil.rmtree(f"{self.settings.workdir}/.scan_temp")
        if self.pymol_animation:
            self.stop_pymol_animation()

        self.react.setup_window = None

    def handle_overall_charge(self, charge):
        # Don't auto-fill charge for xyz files - user must set it manually
        if self.is_xyz_file:
            return

        if self.charge:
            return

        if charge:
            self.charge = charge
            self.ui.lineEdit_charge.setText(self.charge)

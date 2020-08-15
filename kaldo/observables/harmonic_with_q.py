
from kaldo.grid import wrap_coordinates
from kaldo.observables.forceconstant import chi
from kaldo.observables.observable import Observable
import numpy as np
from opt_einsum import contract
import ase.units as units
from kaldo.helpers.storage import lazy_property
import tensorflow as tf
from scipy.linalg.lapack import zheev

from kaldo.helpers.logger import get_logger, log_size
logging = get_logger()

MIN_N_MODES_TO_STORE = 1000
EVTOTENJOVERMOL = units.mol / (10 * units.J)


class HarmonicWithQ(Observable):

    def __init__(self, q_point, second,
                 distance_threshold=None, storage='numpy', is_nw=False, *kargs, **kwargs):
        super().__init__(*kargs, **kwargs)
        self.q_point = q_point
        self.atoms = second.atoms
        self.n_modes = self.atoms.positions.shape[0] * 3
        self.supercell = second.supercell
        self.is_amorphous = (np.array(self.supercell) == [1, 1, 1]).all()
        self.second = second
        self.distance_threshold = distance_threshold
        self.physical_mode= np.ones((1, self.n_modes), dtype=bool)
        self.is_nw = is_nw
        if (q_point == [0, 0, 0]).all():
            if self.is_nw:
                self.physical_mode[0, :4] = False
            else:
                self.physical_mode[0, :3] = False
        if self.n_modes > MIN_N_MODES_TO_STORE:
            self.storage = storage
        else:
            self.storage = 'memory'


    @lazy_property(label='<q_point>')
    def frequency(self):
        frequency = self.calculate_frequency()[np.newaxis, :]
        return frequency


    @lazy_property(label='<q_point>')
    def velocity(self):
        velocity = self.calculate_velocity()
        return velocity


    @lazy_property(label='<q_point>')
    def _dynmat_derivatives_x(self):
        _dynmat_derivatives = self.calculate_dynmat_derivatives_lapack(direction=0)
        return _dynmat_derivatives

    @lazy_property(label='<q_point>')
    def _dynmat_derivatives_y(self):
        _dynmat_derivatives = self.calculate_dynmat_derivatives_lapack(direction=1)
        return _dynmat_derivatives

    @lazy_property(label='<q_point>')
    def _dynmat_derivatives_z(self):
        _dynmat_derivatives = self.calculate_dynmat_derivatives_lapack(direction=2)
        return _dynmat_derivatives


    @lazy_property(label='')
    def _dynmat(self):
        _dynmat = self.calculate_dynmat()
        return _dynmat
    

    @lazy_property(label='<q_point>')
    def _dynmat_fourier(self):
        dynmat_fourier = self.calculate_dynmat_fourier()
        return dynmat_fourier


    @lazy_property(label='<q_point>')
    def _eigensystem(self):
        _eigensystem = self.calculate_eigensystem_lapack(only_eigenvals=False)
        return _eigensystem


    @lazy_property(label='<q_point>')
    def _sij_x(self):
        _sij = self.calculate_sij(direction=0)
        return _sij


    @lazy_property(label='<q_point>')
    def _sij_y(self):
        _sij = self.calculate_sij(direction=1)
        return _sij


    @lazy_property(label='<q_point>')
    def _sij_z(self):
        _sij = self.calculate_sij(direction=2)
        return _sij


    def calculate_frequency(self):
        #TODO: replace calculate_eigensystem() with eigensystem
        eigenvals = self.calculate_eigensystem_lapack(only_eigenvals=True)
        frequency = np.abs(eigenvals) ** .5 * np.sign(eigenvals) / (np.pi * 2.)
        return frequency.real


    def calculate_dynmat_derivatives(self, direction):
        q_point = self.q_point
        is_amorphous = self.is_amorphous
        distance_threshold = self.distance_threshold
        atoms = self.atoms
        list_of_replicas = self.second.list_of_replicas
        replicated_cell = self.second.replicated_atoms.cell
        replicated_cell_inv = self.second._replicated_cell_inv
        cell_inv = self.second.cell_inv
        dynmat = self._dynmat
        positions = self.atoms.positions
        n_unit_cell = atoms.positions.shape[0]
        n_modes = n_unit_cell * 3
        n_replicas = np.prod(self.supercell)
        shape = (1, n_unit_cell * 3, n_unit_cell * 3)
        if is_amorphous:
            type = np.float
        else:
            type = np.complex
        dir = ['_x', '_y', '_z']
        log_size(shape, type, name='dynamical_matrix_derivative_' + dir[direction])
        if is_amorphous:
            distance = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
            distance = wrap_coordinates(distance, replicated_cell, replicated_cell_inv)
            dynmat_derivatives = contract('ij,ibjc->ibjc',
                                          tf.convert_to_tensor(distance[..., direction]),
                                          dynmat[0, :, :, 0, :, :],
                                          backend='tensorflow')
        else:
            distance = positions[:, np.newaxis, np.newaxis, :] - (
                    positions[np.newaxis, np.newaxis, :, :] + list_of_replicas[np.newaxis, :, np.newaxis, :])

            if distance_threshold is not None:

                distance_to_wrap = positions[:, np.newaxis, np.newaxis, :] - (
                    self.second.replicated_atoms.positions.reshape(n_replicas, n_unit_cell, 3)[
                    np.newaxis, :, :, :])

                shape = (n_unit_cell, 3, n_unit_cell, 3)
                type = np.complex
                dynmat_derivatives = np.zeros(shape, dtype=type)
                for l in range(n_replicas):
                    wrapped_distance = wrap_coordinates(distance_to_wrap[:, l, :, :], replicated_cell,
                                                        replicated_cell_inv)
                    mask = (np.linalg.norm(wrapped_distance, axis=-1) < distance_threshold)
                    id_i, id_j = np.argwhere(mask).T
                    dynmat_derivatives[id_i, :, id_j, :, :] += contract('f,fbc->fbc', distance[id_i, l, id_j, direction], \
                                                                         dynmat.numpy()[0, id_i, :, 0, id_j, :] *
                                                                         chi(q_point, list_of_replicas, cell_inv)[l])
            else:

                dynmat_derivatives = contract('ilj,ibljc,l->ibjc',
                                              tf.convert_to_tensor(distance.astype(np.complex)[..., direction]),
                                              tf.cast(dynmat[0], tf.complex128),
                                              tf.convert_to_tensor(chi(q_point, list_of_replicas, cell_inv).flatten().astype(np.complex)),
                                              backend='tensorflow')
        dynmat_derivatives = tf.reshape(dynmat_derivatives, (n_modes, n_modes))
        return dynmat_derivatives


    def calculate_sij(self, direction):
        q_point = self.q_point
        is_amorphous = self.is_amorphous
        shape = (3 * self.atoms.positions.shape[0], 3 * self.atoms.positions.shape[0])
        if is_amorphous:
            type = np.float
        else:
            type = np.complex
        eigenvects = self._eigensystem[1:, :]
        if direction == 0:
            dynmat_derivatives = self._dynmat_derivatives_x
        if direction == 1:
            dynmat_derivatives = self._dynmat_derivatives_y
        if direction == 2:
            dynmat_derivatives = self._dynmat_derivatives_z
        if self.atoms.positions.shape[0] > 100:
            # We want to print only for big systems
            logging.info('Flux operators for q = ' + str(q_point) + ', direction = ' + str(direction))
            dir = ['_x', '_y', '_z']
            log_size(shape, type, name='sij' + dir[direction])
        if is_amorphous:
            sij = tf.tensordot(eigenvects, dynmat_derivatives, (0, 1))
            sij = tf.tensordot(eigenvects, sij, (0, 1))
        else:
            eigenvects = tf.cast(eigenvects, tf.complex128)
            sij = tf.tensordot(eigenvects, dynmat_derivatives, (0, 1))
            sij = tf.tensordot(tf.math.conj(eigenvects), sij, (0, 1))
        return sij


    def calculate_velocity(self):
        frequency = self.frequency[0]
        velocity = np.zeros((self.n_modes, 3))
        inverse_sqrt_freq = tf.cast(tf.convert_to_tensor(1 / np.sqrt(frequency)), tf.complex128)
        for alpha in range(3):
            if alpha == 0:
                sij = self._sij_x
            if alpha == 1:
                sij = self._sij_y
            if alpha == 2:
                sij = self._sij_z
            velocity_AF = 1 / (2 * np.pi) * contract('mn,m,n->mn', sij,
                                   inverse_sqrt_freq, inverse_sqrt_freq, backend='tensorflow') / 2
            velocity_AF = tf.where(tf.math.is_nan(tf.math.real(velocity_AF)), 0., velocity_AF)
            velocity[..., alpha] = contract('mm->m', velocity_AF.numpy().imag)
        return velocity[np.newaxis, ...]


    def calculate_dynmat_fourier(self):
        q_point = self.q_point
        distance_threshold = self.distance_threshold
        atoms = self.atoms
        n_unit_cell = atoms.positions.shape[0]
        n_replicas = np.prod(self.supercell)
        dynmat = self._dynmat
        cell_inv = self.second.cell_inv
        replicated_cell_inv = self.second._replicated_cell_inv
        is_at_gamma = (q_point == (0, 0, 0)).all()
        is_amorphous = (n_replicas == 1)
        list_of_replicas = self.second.list_of_replicas
        log_size((self.n_modes, self.n_modes), np.complex, name='dynmat_fourier')
        if distance_threshold is not None:
            shape = (n_unit_cell, 3, n_unit_cell, 3)
            type = np.complex
            dyn_s = np.zeros(shape, dtype=type)
            replicated_cell = self.second.replicated_atoms.cell

            for l in range(n_replicas):
                distance_to_wrap = atoms.positions[:, np.newaxis, :] - (
                    self.second.replicated_atoms.positions.reshape(n_replicas, n_unit_cell ,3)[np.newaxis, l, :, :])

                distance_to_wrap = wrap_coordinates(distance_to_wrap, replicated_cell, replicated_cell_inv)

                mask = np.linalg.norm(distance_to_wrap, axis=-1) < distance_threshold
                id_i, id_j = np.argwhere(mask).T
                dyn_s[id_i, :, id_j, :] += dynmat.numpy()[0, id_i, :, 0, id_j, :] * chi(q_point, list_of_replicas, cell_inv)[l]
        else:
            if is_at_gamma:
                if is_amorphous:
                    dyn_s = dynmat[0]
                else:
                    dyn_s = contract('ialjb->iajb', dynmat[0], backend='tensorflow')
            else:
                dyn_s = contract('ialjb,l->iajb',
                                 tf.cast(dynmat[0], tf.complex128),
                                 tf.convert_to_tensor(chi(q_point, list_of_replicas, cell_inv).flatten()),
                                 backend='tensorflow')
        dyn_s = tf.reshape(dyn_s, (self.n_modes, self.n_modes))
        return dyn_s


    def calculate_eigensystem(self, only_eigenvals):
        dyn_s = self._dynmat_fourier
        if only_eigenvals:
            esystem = tf.linalg.eigvalsh(dyn_s)
        else:
            log_size(self._dynmat_fourier.shape, type=np.complex, name='eigensystem')
            esystem = tf.linalg.eigh(dyn_s)
            esystem = tf.concat(axis=0, values=(esystem[0][tf.newaxis, :], esystem[1]))
        return esystem


    def calculate_dynmat(self):
        mass = self.atoms.get_masses()
        shape = self.second.value.shape
        log_size(shape, np.float, name='dynmat')
        dynmat = self.second.value * 1 / np.sqrt(mass[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])
        dynmat = dynmat * 1 / np.sqrt(mass[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :, np.newaxis])
        evtotenjovermol = units.mol / (10 * units.J)
        return tf.convert_to_tensor(dynmat * evtotenjovermol)


    def calculate_eigensystem_lapack(self, only_eigenvals=False):
        q_point = self.q_point
        scell = self.supercell
        n_replicas = np.prod(scell)
        atoms = self.atoms
        cell = atoms.cell
        n_unit_cell = atoms.positions.shape[0]
        distance = np.zeros((n_unit_cell, n_unit_cell, 3))
        positions = atoms.positions

        # ev_s = (units._hplanck) * units.J
        # toTHz = 2 * np.pi * units.Rydberg / ev_s * 1e-12
        # massfactor = 2 * units._me * units._Nav * 1000

        # fc_s = forceconstants.second_order.dynmat(atoms.get_masses()) / (Rydberg / (Bohr ** 2))
        # EVTOTENJOVERMOL = units.mol / (10 * units.J)
        # fc_s = fc_s / EVTOTENJOVERMOL * massfactor
        fc_s = self._dynmat.numpy()
        replicated_positions = self.second.replicated_atoms.positions.reshape((n_replicas, n_unit_cell, 3))

        # mass = self.atoms.get_masses()
        shape = fc_s.shape
        log_size(shape, np.float, name='dynmat')
        # fc_s = second * 1 / np.sqrt(mass[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis, np.newaxis])
        # fc_s = fc_s * 1 / np.sqrt(mass[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :, np.newaxis])
        # evtotenjovermol = units.mol / (10 * units.J)
        # fc_s = (fc_s * evtotenjovermol)
        fc_s = fc_s.reshape((n_unit_cell, 3, scell[0], scell[1], scell[2], n_unit_cell, 3))

        replicated_cell = cell * scell
        ir = 0


        sc_r_pos = np.zeros((3 ** 3, 3))
        sc_r_pos_norm = np.zeros((3 ** 3))

        for ix2 in [-1, 0, 1]:
            for iy2 in [-1, 0, 1]:
                for iz2 in [-1, 0, 1]:

                    for i in np.arange(3):
                        sc_r_pos[ir, i] = np.dot(replicated_cell[:, i], np.array([ix2,iy2,iz2]))
                    sc_r_pos_norm[ir] = 0.5 * np.dot(sc_r_pos[ir], sc_r_pos[ir])
                    ir = ir + 1
        dyn_s = np.zeros((n_unit_cell, 3, n_unit_cell, 3), dtype=np.complex)
        ddyn_s = np.zeros((n_unit_cell, 3, n_unit_cell, 3, 3), dtype=np.complex)
        list_of_index = np.round((replicated_positions - self.atoms.positions).dot(
            np.linalg.inv(atoms.cell))).astype(np.int)
        list_of_index = list_of_index[:, 0, :]

        tt = []
        rreplica = []
        for ix2 in [-1, 0, 1]:
            for iy2 in [-1, 0, 1]:
                for iz2 in [-1, 0, 1]:
                    for f in range(list_of_index.shape[0]):

                        scell_id = np.array([ix2 * scell[0], iy2 * scell[1], iz2 * scell[2]])
                        replica_id = list_of_index[f]
                        t = replica_id + scell_id
                        replica_position = np.tensordot(cell, t, (0, -1))
                        tt.append(t)
                        rreplica.append(replica_position)

        tt = np.array(tt)
        rreplica = np.array(rreplica)
        for ind in range(tt.shape[0]):
            t = tt[ind]
            replica_position = rreplica[ind]
            for iat in np.arange(n_unit_cell):
                for jat in np.arange(n_unit_cell):
                    distance = replica_position + (positions[iat, :] - positions[jat, :])
                    projection = (np.dot(sc_r_pos, distance) - sc_r_pos_norm[:])

                    if ((projection <= 1e-6).all()):
                        neq = (np.abs(projection) <= 1e-6).sum()
                        weight = 1.0 / (neq)
                        qr = 2. * np.pi * np.dot(q_point[:], t[:])
                        for ipol in np.arange(3):
                            for jpol in np.arange(3):
                                dyn_s[iat, ipol, jat, jpol] += fc_s[
                                     jat, jpol, t[0], t[1], t[2], iat, ipol] * np.exp(-1j * qr) * weight


        frequency = np.zeros((n_unit_cell * 3))
        if only_eigenvals:
            esystem = np.zeros((n_unit_cell * 3), dtype=np.complex)
        else:
            esystem = np.zeros((n_unit_cell * 3 + 1, n_unit_cell * 3), dtype=np.complex)

        dyn = dyn_s[...].reshape((n_unit_cell * 3, n_unit_cell * 3))

        omega2,eigenvect,info = zheev(dyn)
        frequency[:] = np.sign(omega2) * np.sqrt(np.abs(omega2))
        frequency[:] = frequency[:] / np.pi / 2
        if only_eigenvals:
            esystem = (frequency[ :] * np.pi * 2) ** 2
        else:
            esystem = np.vstack(((frequency[:] * np.pi * 2) ** 2, eigenvect))
        return esystem


    def calculate_dynmat_derivatives_lapack(self, direction=None):
        # Debugging method to compare results, do not use
        q_point = self.q_point
        supercell = self.supercell
        atoms = self.atoms
        cell = atoms.cell
        n_unit_cell = atoms.positions.shape[0]
        replicated_positions = self.second.replicated_atoms.positions
        # mass = atoms.get_masses()
        # masses_2d = np.zeros((n_unit_cell, n_unit_cell))
        distance = np.zeros((n_unit_cell, n_unit_cell, 3))
        positions = atoms.positions
        replicated_cell = cell * supercell
        fc_s = self._dynmat.numpy()
        # fc_s = fc_s / (Rydberg / (Bohr ** 2))
        # EVTOTENJOVERMOL = units.mol / (10 * units.J)
        # massfactor = 2 * units._me * units._Nav * 1000
        # fc_s = fc_s / EVTOTENJOVERMOL * massfactor
        fc_s = fc_s.reshape((n_unit_cell, 3, supercell[0], supercell[1], supercell[2], n_unit_cell, 3))
        j = 0

        sc_r_pos = np.zeros((3 ** 3, 3))
        ir = 0
        for ix2 in [-1, 0, 1]:
            for iy2 in [-1, 0, 1]:
                for iz2 in [-1, 0, 1]:

                    for i in np.arange(3):
                        sc_r_pos[ir, i] = np.dot(replicated_cell[:, i], np.array([ix2,iy2,iz2]))
                    ir = ir + 1

        sc_r_pos_norm = 1 / 2 * np.linalg.norm(sc_r_pos, axis=1) ** 2
        ddyn_s = np.zeros((n_unit_cell, 3, n_unit_cell, 3, 3), dtype=np.complex)
        list_of_index = np.round((replicated_positions - atoms.positions).dot(
            np.linalg.inv(atoms.cell))).astype(np.int)
        list_of_index = list_of_index[:, 0, :]
        tt = []
        rreplica = []
        for ix2 in [-1, 0, 1]:
            for iy2 in [-1, 0, 1]:
                for iz2 in [-1, 0, 1]:
                    for f in range(list_of_index.shape[0]):

                        scell_id = np.array([ix2 * supercell[0], iy2 * supercell[1], iz2 * supercell[2]])
                        replica_id = list_of_index[f]
                        t = replica_id + scell_id
                        replica_position = np.tensordot(cell, t, (0, -1))
                        tt.append(t)
                        rreplica.append(replica_position)

        tt = np.array(tt)
        rreplica = np.array(rreplica)
        for ind in range(tt.shape[0]):
            t = tt[ind]
            replica_position = rreplica[ind]

            for iat in np.arange(n_unit_cell):
                for jat in np.arange(n_unit_cell):
                    distance = replica_position + (positions[iat] - positions[jat])
                    projection = (np.dot(sc_r_pos, distance) - sc_r_pos_norm)

                    if (projection <= 1e-6).all():
                        neq = (np.abs(projection) <= 1e-6).sum()
                        weight = 1.0 / (neq)
                        qr = 2. * np.pi * np.dot(q_point[:], t[:])
                        for ipol in np.arange(3):
                            for jpol in np.arange(3):
                                ddyn_s[iat, ipol, jat, jpol, :] -= replica_position * fc_s[
                                    jat, jpol, t[0], t[1], t[2], iat, ipol] * np.exp(-1j * qr) * weight
        return ddyn_s.reshape((n_unit_cell * 3, n_unit_cell * 3, 3))



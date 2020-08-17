
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
                 distance_threshold=None, storage='numpy', is_nw=False, is_unfolding=False, *kargs, **kwargs):
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
        self.is_unfolding = is_unfolding
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
        if self.is_unfolding:
            _dynmat_derivatives = self.calculate_dynmat_derivatives_unfolded(direction=0)
        else:
            _dynmat_derivatives = self.calculate_dynmat_derivatives(direction=0)
        return _dynmat_derivatives


    @lazy_property(label='<q_point>')
    def _dynmat_derivatives_y(self):
        if self.is_unfolding:
            _dynmat_derivatives = self.calculate_dynmat_derivatives_unfolded(direction=1)
        else:
            _dynmat_derivatives = self.calculate_dynmat_derivatives(direction=1)
        return _dynmat_derivatives


    @lazy_property(label='<q_point>')
    def _dynmat_derivatives_z(self):
        if self.is_unfolding:
            _dynmat_derivatives = self.calculate_dynmat_derivatives_unfolded(direction=2)
        else:
            _dynmat_derivatives = self.calculate_dynmat_derivatives(direction=2)
        return _dynmat_derivatives
    

    @lazy_property(label='<q_point>')
    def _dynmat_fourier(self):
        dynmat_fourier = self.calculate_dynmat_fourier()
        return dynmat_fourier


    @lazy_property(label='<q_point>')
    def _eigensystem(self):
        if self.is_unfolding:
            _eigensystem = self.calculate_eigensystem_unfolded(only_eigenvals=False)
        else:
            _eigensystem = self.calculate_eigensystem(only_eigenvals=False)
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
        if self.is_unfolding:
            eigenvals = self.calculate_eigensystem_unfolded(only_eigenvals=True)
        else:
            eigenvals = self.calculate_eigensystem(only_eigenvals=True)
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
        dynmat = self.second.dynmat
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
        dynmat = self.second.dynmat
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


    def calculate_eigensystem_unfolded(self, only_eigenvals=False):
        q_point = self.q_point
        scell = self.supercell
        atoms = self.atoms
        cell = atoms.cell
        n_unit_cell = atoms.positions.shape[0]
        positions = atoms.positions
        fc_s = self.second.dynmat.numpy()
        fc_s = fc_s.reshape((n_unit_cell, 3, scell[0], scell[1], scell[2], n_unit_cell, 3))
        sc_r_pos = self.second.supercell_positions
        sc_r_pos_norm = 1 / 2 * np.linalg.norm(sc_r_pos, axis=1) ** 2
        dyn_s = np.zeros((n_unit_cell, 3, n_unit_cell, 3), dtype=np.complex)
        tt = self.second.supercell_replicas
        for ind in range(tt.shape[0]):
            t = tt[ind]
            replica_position = np.tensordot(t, cell, (-1, 0))
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
        dyn = dyn_s[...].reshape((n_unit_cell * 3, n_unit_cell * 3))
        omega2,eigenvect,info = zheev(dyn)
        frequency = np.sign(omega2) * np.sqrt(np.abs(omega2))
        frequency = frequency[:] / np.pi / 2
        if only_eigenvals:
            esystem = (frequency[:] * np.pi * 2) ** 2
        else:
            esystem = np.vstack(((frequency[:] * np.pi * 2) ** 2, eigenvect))
        return esystem


    def calculate_dynmat_derivatives_unfolded(self, direction=None):
        q_point = self.q_point
        supercell = self.supercell
        atoms = self.atoms
        cell = atoms.cell
        n_unit_cell = atoms.positions.shape[0]
        ddyn_s = np.zeros((n_unit_cell, 3, n_unit_cell, 3), dtype=np.complex)
        positions = atoms.positions
        fc_s = self.second.dynmat.numpy()
        fc_s = fc_s.reshape((n_unit_cell, 3, supercell[0], supercell[1], supercell[2], n_unit_cell, 3))
        sc_r_pos = self.second.supercell_positions
        sc_r_pos_norm = 1 / 2 * np.linalg.norm(sc_r_pos, axis=1) ** 2
        tt = self.second.supercell_replicas
        for ind in range(tt.shape[0]):
            t = tt[ind]
            replica_position = np.tensordot(t, cell, (-1, 0))
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
                                ddyn_s[iat, ipol, jat, jpol] -= replica_position[direction] * fc_s[
                                    jat, jpol, t[0], t[1], t[2], iat, ipol] * np.exp(-1j * qr) * weight
        return ddyn_s.reshape((n_unit_cell * 3, n_unit_cell * 3))



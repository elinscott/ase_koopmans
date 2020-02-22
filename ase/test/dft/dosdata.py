import matplotlib.pyplot as plt
import numpy as np
import pytest

from ase.dft.dosdata import DOSData, GridDOSData, RawDOSData


class MinimalDOSData(DOSData):
    """Inherit from ABC to test its features"""
    def get_energies(self):
        super().get_energies()

    def get_weights(self):
        super().get_weights()


class TestDosData:
    """Test the abstract base class for DOS data"""

    sample_info = [(None, {}),
                   ({}, {}),
                   ({'symbol': 'C', 'index': '2', 'strangekey': 'isallowed'},
                    {'symbol': 'C', 'index': '2', 'strangekey': 'isallowed'}),
                   ('notadict', TypeError),
                   (False, TypeError)]

    @pytest.mark.parametrize('info, expected', sample_info)
    def test_dosdata_init_info(self, info, expected):
        """Check 'info' parameter is handled properly"""
        if isinstance(expected, type) and isinstance(expected(), Exception):
            with pytest.raises(expected):
                dos_data = MinimalDOSData(info=info)
        else:
            dos_data = MinimalDOSData(info=info)
            assert dos_data.info == expected

    dosdata_abc_notimplemented_methods_args = [('get_energies', tuple()),
                                               ('get_weights', tuple())]
    @pytest.mark.parametrize('method, args',
                             dosdata_abc_notimplemented_methods_args)
    def test_dosdata_notimplemented(self, method, args):
        """Check NotImplementedError raised from abstract base class"""
        dos_data = MinimalDOSData()
        with pytest.raises(NotImplementedError):
            getattr(dos_data, method)(*args)


class TestRawDosData:
    """Test the raw DOS data container"""

    @pytest.fixture
    def sparse_dos(self):
        return RawDOSData([1.2, 3.4, 5.], [3., 2.1, 0.],
                          info={'symbol': 'H', 'number': '1', 'food': 'egg'})

    @pytest.fixture
    def another_sparse_dos(self):
        return RawDOSData([8., 2., 2., 5.], [1., 1., 1., 1.],
                          info={'symbol': 'H', 'number': '2'})

    def test_access(self, sparse_dos):
        assert sparse_dos.info == {'symbol': 'H', 'number': '1', 'food': 'egg'}
        assert np.allclose(sparse_dos.get_energies(), [1.2, 3.4, 5.])
        assert np.allclose(sparse_dos.get_weights(), [3., 2.1, 0.])

    def test_addition(self, sparse_dos, another_sparse_dos):
        summed_dos = sparse_dos + another_sparse_dos
        assert summed_dos.info == {'symbol': 'H'}
        assert np.allclose(summed_dos.get_energies(),
                           [1.2, 3.4, 5., 8., 2., 2., 5.])
        assert np.allclose(summed_dos.get_weights(),
                           [3., 2.1, 0., 1., 1., 1., 1.])

    sampling_data_args_results = [
        # Special case: peak max at width 1
        ([[0.], [1.]],
         [[0.], {'width': 1}],
         [1. / (np.sqrt(2. * np.pi))]),
        # Peak max with different width, position
        ([[1.], [2.]],
         [[1.], {'width': 0.5}],
         [2. / (np.sqrt(2. * np.pi) * 0.5)]),
        # Peak max for two simultaneous deltas
        ([[1., 1.], [2., 1.]],
         [[1.], {'width': 1}],
         [3. / (np.sqrt(2. * np.pi))]),
        # Compare with theoretical half-maximum
        ([[0.], [1.]],
         [[np.sqrt(2 * np.log(2)) * 3],
          {'width': 3}],
         [0.5 / (np.sqrt(2 * np.pi) * 3)]),
        # And a case with multiple values, generated
        # using the ASE code (not benchmarked)
        ([[1.2, 3.4, 5], [3., 2.1, 0.]],
         [[1., 1.5, 2., 2.4], {'width': 2}],
         [0.79932418, 0.85848101, 0.88027184, 0.8695055])]

    @pytest.mark.parametrize('data, args, result',
                             sampling_data_args_results)
    def test_sampling(self, data, args, result):
        dos = RawDOSData(data[0], data[1])
        assert np.allclose(dos.sample(*args[:-1], **args[-1]), result)
        
    def test_sampling_error(self, sparse_dos):
        with pytest.raises(ValueError):
            sparse_dos.sample([1, 2, 3], width=0.)
        with pytest.raises(ValueError):
            sparse_dos.sample([1, 2, 3], width=-1)

    def test_sample_grid(self, sparse_dos):
        min_dos = sparse_dos.sample_grid(10, xmax=5, padding=3, width=0.1)
        assert min_dos[0][0] == 1.2 - 3 * 0.1

        max_dos = sparse_dos.sample_grid(10, xmin=0, padding=2, width=0.2)
        assert max_dos[0][-1] == 5 + 2 * 0.2

        default_dos = sparse_dos.sample_grid(10)
        assert np.allclose(default_dos[0], np.linspace(0.9, 5.3, 10))
        assert np.allclose(default_dos[1],
                           sparse_dos.sample(np.linspace(0.9, 5.3, 10)))

    # Comparing plot outputs is hard, so we
    # - inspect the line values
    # - check that a line styling parameter is correctly passed through mplargs
    # - set a kwarg from self.sample() to check broadening args are recognised
    linewidths = [1, 5]
    @pytest.mark.parametrize('linewidth', linewidths)
    def test_plot_dos(self, sparse_dos, linewidth):
        fig, ax = plt.subplots()
        sparse_dos.plot_dos(npts=5, ax=ax,
                            mplargs={'linewidth': linewidth}, smearing='Gauss')

        assert ax.lines[0].get_linewidth() == linewidth
        line_data = ax.lines[0].get_data()
        assert np.allclose(line_data[0], np.linspace(0.9, 5.3, 5))
        assert np.allclose(line_data[1],
                           [1.32955452e-01, 1.51568133e-13,
                            9.30688167e-02, 1.06097693e-13, 3.41173568e-78])

    @pytest.mark.parametrize('linewidth', linewidths)
    def test_plot_deltas(self, sparse_dos, linewidth):
        fig, ax = plt.subplots()
        sparse_dos.plot_deltas(ax=ax, mplargs={'linewidth': linewidth})

        assert ax.get_children()[0].get_linewidth() == linewidth

        assert np.allclose(list(map(lambda x: x.vertices,
                                    ax.get_children()[0].get_paths())),
                           [[[1.2, 0.], [1.2, 3.]],
                            [[3.4, 0.], [3.4, 2.1]],
                            [[5., 0.], [5., 0.]]])


class TestGridDosData:
    """Test the grid DOS data container"""
    def test_init(self):
        # energies and weights must be equal lengths
        with pytest.raises(ValueError):
            GridDOSData(np.linspace(0, 10, 11), np.zeros(10))

        # energies must be evenly spaced
        with pytest.raises(ValueError):
            GridDOSData(np.linspace(0, 10, 11)**2, np.zeros(11))

    @pytest.fixture
    def dense_dos(self):
        x = np.linspace(0., 10., 11)
        y = np.sin(x / 10)
        return GridDOSData(x, y, info={'symbol': 'C', 'orbital': '2s',
                                       'day': 'Tue'})

    @pytest.fixture
    def another_dense_dos(self):
        x = np.linspace(0., 10., 11)
        y = np.sin(x / 10) * 2
        return GridDOSData(x, y, info={'symbol': 'C', 'orbital': '2p',
                                       'month': 'Feb'})

    def test_access(self, dense_dos):
        assert dense_dos.info == {'symbol': 'C', 'orbital': '2s', 'day': 'Tue'}
        assert len(dense_dos.get_energies()) == 11
        assert dense_dos.get_energies()[-2] == 9.
        assert dense_dos.get_weights()[-1] == np.sin(1)

    def test_addition(self, dense_dos, another_dense_dos):
        sum_dos = dense_dos + another_dense_dos
        assert np.allclose(sum_dos.get_energies(), dense_dos.get_energies())
        assert np.allclose(sum_dos.get_weights(), dense_dos.get_weights() * 3)
        assert sum_dos.info == {'symbol': 'C'}

        with pytest.raises(TypeError):
            dense_dos + RawDOSData([1., 2.], [3., 4.])

        with pytest.raises(ValueError):
            dense_dos + GridDOSData(dense_dos.get_energies()[1:],
                                    dense_dos.get_weights()[1:])

    @pytest.mark.parametrize('npts, width, padding', [(19, 1, 4)])
    def test_grid_sampling(self, dense_dos, npts, width, padding):
        x, y = dense_dos.sample_grid(npts, width=width, padding=padding)

        # Check padding for auto limits
        assert x[0] == 0 - width * padding
        assert x[-1] == 10 + width * padding

        assert np.allclose(y, dense_dos.sample(np.linspace(x[0], x[-1], npts),
                                               width=width), atol=1e-3, rtol=10e-2)

import numpy


class PLAMSMatrix(numpy.ndarray):
    """
    PLAMS alternative to a sympy matrix

    Note: Here and there I round decimals to 1.e-10. Marked with comment 'Rounding'
    """

    def __new__(cls, *args, **kwargs):
        """
        Creates the instance
        """
        matrix = args[0]
        if not isinstance(matrix, numpy.ndarray):
            raise Exception("ndarray or ndarray subclass required as argument.")
        if len(args) > 1:
            raise Exception("PLAMSMatrix requires a numpy.ndarray as the single argument.")
        elif len(kwargs) > 0:
            raise Exception("PLAMSMatrix requires a numpy.ndarray as the single argument.")
        return matrix.view(cls)

    def __array_finalize__(self, obj):
        """
        Now initialize the instance
        """
        # FIXME: Slicing should not return a PLAMSMatrix object if the dimension of the slice is 1
        if len(self.shape) > 2:
            raise Exception("A PLAMSMatrix can only be of dimension 2 ", self.shape)

    def nullspace(self):
        """
        Get the nullspace vectors using Gaussian elimination
        """
        n = self.shape[0]
        rank = numpy.linalg.matrix_rank(self)
        # print ('Rank: ',rank)
        # print ('NElements: ',n)

        # Create the big matrix
        bigmat = self._get_big_matrix()
        matrix = bigmat.transpose()
        # print ('Matrix to be eliminated')
        # print(matrix.as_string())

        matrix = matrix.rref(rank, n)
        # print ('Row echolon matrix')
        # print(matrix.as_string())

        # Get the basis vectors
        ut = matrix.transpose()
        # print ('Transposed matrix')
        # print (ut.as_string())
        basis = ut[n:, rank:].transpose()

        # There seems to be some numerical noise that constitutes the difference with the sympy result
        # Not sure why. Rounding to arbitrary 10 decimals
        basis = basis.round(decimals=10)

        return basis

    def rref(self, max_row=None, max_col=None):
        """
        Reduce self to row-echolon form

        * ``max_row`` -- The row up to which to diagonalize
        * ``max_col``    -- The column up to which to put the rows below max_row to zero
        """
        # Set max_row and max_col to appropriate values, if unset
        rank = numpy.linalg.matrix_rank(self)
        ncols = self.shape[1]
        if max_row is None:
            max_row = rank
        if max_col is None:
            max_col = ncols

        # Create the new matrix
        matrix = self.copy()

        # Now we can start doing row elimination
        matrix._perform_gaussian_elimination(max_row, max_col)

        # At this point we should have successfully eliminated up to the max_row of the matrix
        # print ("Diagonalized matrix")
        # print(matrix.as_string())

        # Now the remaining columns can be fully set to zero, at least up to max_col
        matrix._gaussian_elimination_of_region(max_row, max_col)

        return matrix

    def get_heading_zeros(self, irow):
        """
        Get the number of zeros at the start of this row
        """
        row = self[irow]
        ncols = row.shape[0]
        nzeros = 0
        for i in range(ncols):
            if row[i] != 0:
                break
            nzeros += 1
        return nzeros

    def as_string(self, space=8):
        """
        Print a numpy matrix in nice format
        """
        form = "%{}s".format(space)
        form2 = "%{}.1e".format(space)
        lines = []
        for row in self:
            strings = [form % (str(v)) if len(str(v)) <= space else form2 % (v) for v in row]
            s = " ".join(strings)
            lines.append(s)
        return "\n".join(lines)

    #################
    # Private methods
    #################

    def _get_big_matrix(self):
        """
        Get the big matrix upon which we want to do column elinimation
        """
        n = self.shape[0]
        m = self.shape[1]
        bigmat = PLAMSMatrix(numpy.zeros((n + m, m)))
        for i in range(m):
            bigmat[n + i, i] = 1
        bigmat[:n] = self
        return bigmat

    def _perform_gaussian_elimination(self, max_row, max_col):
        """
        Perform Gaussian elimination up to max_row
        """
        nrows = self.shape[0]
        for irow in range(max_row):
            # If there is no suitable row the first time around,
            # try to find one by shifting icol to the right.
            found = False
            for icol in range(irow, max_col):
                for iswap in range(nrows + 1 - irow):
                    A = self.copy()
                    A._set_row_to_echolon_form(irow, icol)
                    # The intension here is that self[irow,icol] becomes zero
                    # Rounding to arbitrary 10 decimals
                    A[irow][abs(A[irow]) < 1e-10] = 0
                    # If there are too many zeros, move the row to the end, and start again
                    nzeros = self.get_heading_zeros(irow)
                    nzerosA = A.get_heading_zeros(irow)
                    if nzeros > irow or nzerosA > irow:
                        self._shift_row_to(irow, nrows - 1)
                    else:
                        found = True
                        break
                if found:
                    break
            # If no good row was found, we are back to the original one, and set that one
            self[:, :] = A

    def _gaussian_elimination_of_region(self, min_row=None, max_col=None):
        """
        Perform Gaussian elimination only on the rows min_row and higher, and only up to max_col

        Note: The region is selected so that all values can be set to zero
        """
        nrows = self.shape[0]
        icount = min_row
        for irow in range(min_row, nrows):
            # For each of these rows, put the first max_col values to zero
            self._set_row_to_echolon_form(icount, max_col)
            # Rounding to arbitray 10 decimals
            self[icount][abs(self[icount]) < 1e-10] = 0
            if (self[icount, :max_col] == 0).all() and icount < nrows - 1:
                self._shift_row_to(icount, nrows - 1)
            else:
                icount += 1

    def _set_row_to_echolon_form(self, irow, final=None):
        """
        Eliminate values in row irow, with zeros up to column final

        FIXME: The final variable is not used. Change?
        """
        if final is None:
            final = irow
        # for icol in range(irow):
        for icol in range(final):
            # Subtract a multiple of the row icol from this row
            prev_col = icol
            if icol >= irow:
                # then check if the previous row has enough zeros
                prev_col = irow - 1
                nzeros = self.get_heading_zeros(prev_col)
                if nzeros < icol:
                    continue
            if self[prev_col, icol] == 0:
                continue
            scale = -self[irow, icol] / self[prev_col, icol]
            self._row_add(prev_col, irow, scale)

    def _shift_row_to(self, k, l):
        """
        Shift row k to a later row l, and have the rest shift up
        """
        if l <= k:
            raise Exception("Bad choice of rows")
        row = self[k].copy()
        self[k:l] = self[k + 1 : l + 1]
        self[l] = row

    def _row_add(self, k, l, scale):
        """
        Changes row l: Adds values of row k multiplied by scale.
        """
        self[l] += self[k] * scale

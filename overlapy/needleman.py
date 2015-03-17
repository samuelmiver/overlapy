import numpy as np

class Needleman:
    """
    This class runs the Needleman & Wunsch for two given sequences (lists) of objects, a score matrix and
    some gap penalty. It returns a tuple with two lists of the same objects that you passed on input or none
    if there is a gap.

    It follows the extract method to method object design pattern. You enter all the input parameters,
    then you call align() and you collect the results. You will have to create a new class
    for each alignment. YOU CAN'T REUSE INSTANCES.

    Sample usage:

        >>> nw = Needleman(seq1, seq2, score_matrix)
        >>> nw.align()
        >>> alignment = nw.get_alignment()
        >>> score = nw.score

    """

    def __init__(self, sequence1, sequence2, score_matrix, gap_penalty = 10):
        """
        Given two lists of items (eg. nucleotides, aminoacids, etc.) and a
        score matrix, perform an alignment of the two lists and return the
        result as a couple of iterables. Each of them will contain the aligned
        original objects, with gaps represented by elements of type None.
        """
        self.sequences = [sequence1, sequence2]
        self.score_matrix = score_matrix
        self.gap_penalty = gap_penalty
        self.already_executed = False

    def align(self):
        """
        Given the correct input has been given, aligns the two sequences and sets the return values.
        """
        if not self.already_executed:
            pointers = self.__get_pointers()
            self.alignment = self.__calculate_alignment(pointers)
            self.already_executed = True

    def get_alignment(self):
        """
        Returns an alignment, which is comprised of two lists in a tuple. If there is a gap it will be
        represented by None.
        """
        return self.alignment

    # PRIVATE

    def __get_pointers(self):
        """
        Given the initial matrix of scores (type np.matrix), and the gap
        penalty, calculate the matrix of accumulated scores and the matrix
        of pointers. Return the matrix of pointers (type np.matrix).
        """

        # we need an extra row and col for the gaps
        rows = self.score_matrix.shape[0] + 1
        cols = self.score_matrix.shape[1] + 1

        accumulated_scores = np.empty([rows, cols])
        # matrix of pointers, with the same dimensions
        # 0 = diagonal pointer, 1 = pointer to top, 2 = pointer to left
        pointers = np.empty([rows, cols])

        # fill first column
        for i in range(rows):
            accumulated_scores[i][0] = i*self.gap_penalty
            pointers[i][0] = 1 # point to top

        # fill first row
        for j in range(cols):
            accumulated_scores[0][j] = j*self.gap_penalty
            pointers[0][j] = 2

        # fill the rest of the matrix
        for i in range(1, rows):
            for j in range(1, cols):
                values = []
                # score to go diagonally (pointer = 0)
                values.append(accumulated_scores[i-1][j-1] + self.score_matrix[i-1][j-1])
                # score to go horizontally (pointer = 1)
                values.append(self.score_matrix[i-2][j-1] + self.gap_penalty)
                # score to go vertically (pointer = 2)
                values.append(self.score_matrix[i-1][j-2] + self.gap_penalty)
                min_value = min(values)

                # fill matrix with the best score and set pointer
                accumulated_scores[i][j] = min_value
                pointers[i][j] = values.index(min_value)

        self.score = accumulated_scores[rows-1][cols-1]
        return pointers

    def __calculate_alignment(self, pointers):
        """
        Given two lists of sequences and a matrix of pointers (obtained
        from __get_pointers()), return the aligned sequences as a couple of
        iterables. Each of them will contain the aligned original objects,
        with gaps represented by elements of type None.

        """

        # traceback
        alignment1 = []
        alignment2 = []

        i = pointers.shape[0] - 1
        j = pointers.shape[1] - 1


        while not (i == 0 and j == 0):
            # move in the diagonal?
            pointer= pointers[i][j]
            if pointer == 0:
                alignment1.insert(0, self.sequences[0][i-1])
                alignment2.insert(0, self.sequences[1][j-1])
                i -= 1
                j -= 1
            # move horizontally
            elif pointer == 1:
                alignment1.insert(0, self.sequences[0][i-1])
                # None = gap
                alignment2.insert(0, None)
                i -= 1
            # move vertically
            elif pointer == 2:
                alignment1.insert(0, None)
                alignment2.insert(0, self.sequences[1][j-1])
                j -= 1
            else:
                raise Exception('Unknown pointer value: {} is not 0, 1 or 2'.format(pointer))

        return alignment1, alignment2


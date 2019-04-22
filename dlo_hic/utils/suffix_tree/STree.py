ASCII_SYMBOLS = ['$', '#'] + [chr(i) for i in range(45, 58)]
UNICODE_SYMBOLS = [chr(i) for i in list(range(0xE000, 0xF8FF + 1)) + list(range(0xF0000, 0xFFFFD + 1)) + list(
    range(0x100000, 0x10FFFD + 1))]


class STree():
    """Class representing the suffix tree."""

    def __init__(self, input='', end_symbol_set='unicode'):
        self.root = _SNode()
        self.root.depth = 0
        self.root.idx = 0
        self.root.parent = self.root
        self.root._add_suffix_link(self.root)

        self.words = None
        self.end_symbol_set = end_symbol_set
        self._check_symbol_set(input)

        if input:
            self.build(input)

    def _check_symbol_set(self, input):
        if self.end_symbol_set == 'ascii':
            for c in ASCII_SYMBOLS:
                if c in input:
                    self.end_symbol_set = 'unicode'
                    break

    def _check_input(self, input):
        """Checks the validity of the input.

        In case of an invalid input throws ValueError.
        """
        if isinstance(input, str):
            return 'st'
        elif isinstance(input, list):
            if all(isinstance(item, str) for item in input):
                return 'gst'

        raise ValueError("String argument should be of type String or"
                         " a list of strings")

    def build(self, x):
        """Builds the Suffix tree on the given input.
        If the input is of type List of Strings:
        Generalized Suffix Tree is built.

        :param x: String or List of Strings
        """
        type_ = self._check_input(x)

        if type_ == 'st':
            x += next(self._terminalSymbolsGenerator())
            self._build(x)
        if type_ == 'gst':
            self._build_generalized(x)

    def _build(self, x, algorithm="McCreight"):
        """Builds a Suffix tree."""
        self.word = x
        if algorithm == "McCreight":
            self._build_McCreight(x)
        elif algorithm == "Ukkonen":
            self._build_Ukkonen(x)
        else:
            raise NotImplementedError()

    def _build_McCreight(self, x):
        """Builds a Suffix tree using McCreight O(n) algorithm.

        Algorithm based on:
        McCreight, Edward M. "A space-economical suffix tree construction algorithm." - ACM, 1976.
        Implementation based on:
        UH CS - 58093 String Processing Algorithms Lecture Notes
        """
        node = self.root
        d = 0
        for i in range(len(x)):
            while node.depth == d and node._has_transition(x[d + i]):
                node = node._get_transition_link(x[d + i])
                d = d + 1
                while d < node.depth and x[node.idx + d] == x[i + d]:
                    d = d + 1
            if d < node.depth:
                node = self._insert_node(x, node, d)
            self._create_leaf(x, i, node, d)
            if not node._get_suffix_link():
                self._compute_slink(x, node)
            node = node._get_suffix_link()
            d = d - 1
            if d < 0:
                d = 0

    def _insert_node(self, input, node, depth):
        """
        Insert node between node and node.parent

        node.parent -> new_node -> node

        """
        i = node.idx
        p = node.parent
        new_node = _SNode(idx=i, depth=depth)
        new_node._add_transition_link(node, input[i + depth])
        node.parent = new_node
        p._add_transition_link(new_node, input[i + p.depth])
        new_node.parent = p
        return new_node

    def _create_leaf(self, input, idx, node, depth):
        """
        create a leaf node

        node -> leaf
        """
        leaf = _SNode()
        leaf.idx = idx
        leaf.depth = len(input) - idx
        edge_label = input[idx + depth]
        node._add_transition_link(leaf, edge_label)
        leaf.parent = node
        return leaf

    def _compute_slink(self, input, node):
        """
        """
        d = node.depth
        v = node.parent._get_suffix_link()
        while v.depth < d - 1:
            v = v._get_transition_link(input[node.idx + v.depth + 1])
        if v.depth > d - 1:
            v = self._insert_node(input, v, d - 1)
        node._add_suffix_link(v)

    def _build_Ukkonen(self, x):
        """Builds a Suffix tree using Ukkonen's online O(n) algorithm.

        Algorithm based on:
        Ukkonen, Esko. "On-line construction of suffix trees." - Algorithmica, 1995.
        """
        raise NotImplementedError()

    def _build_generalized(self, xs):
        """Builds a Generalized Suffix Tree (GST) from the array of strings provided.
        """
        terminal_gen = self._terminalSymbolsGenerator()

        _xs = [x + next(terminal_gen) for x in xs]
        self.words = _xs
        _xs = ''.join(_xs)
        self.word = _xs
        self._generalized_word_starts(xs)
        self._build(_xs)
        self.root._traverse(self._label_generalized)

    def _label_generalized(self, node):
        """Helper method that labels the nodes of GST with indexes of strings
        found in their descendants.
        """
        if node.is_leaf():
            x = {self._get_word_start_index(node.idx)}
        else:
            x = {n for ns in node.transition_links for n in ns[0].generalized_idxs}
        node.generalized_idxs = x

    def _get_word_start_index(self, idx):
        """Helper method that returns the index of the string based on node's
        starting index"""
        i = 0
        for _idx in self.word_starts[1:]:
            if idx < _idx:
                return i
            else:
                i += 1
        return i

    def lcs(self, stringIdxs=-1):
        """Returns the Largest Common Substring of Strings provided in stringIdxs.
        If stringIdxs is not provided, the LCS of all strings is returned.

        ::param stringIdxs: Optional: List of indexes of strings.
        """
        if stringIdxs == -1 or not isinstance(stringIdxs, list):
            stringIdxs = set(range(len(self.word_starts)))
        else:
            stringIdxs = set(stringIdxs)

        deepestNode = self._find_lcs(self.root, stringIdxs)
        start = deepestNode.idx
        end = deepestNode.idx + deepestNode.depth
        return self.word[start:end]

    def _find_lcs(self, node, stringIdxs):
        """Helper method that finds LCS by traversing the labeled GSD."""
        nodes = [self._find_lcs(n, stringIdxs)
                 for (n, _) in node.transition_links
                 if n.generalized_idxs.issuperset(stringIdxs)]

        if nodes == []:
            return node

        deepestNode = max(nodes, key=lambda n: n.depth)
        return deepestNode

    def _generalized_word_starts(self, xs):
        """Helper method returns the starting indexes of strings in GST"""
        self.word_starts = []
        i = 0
        for n in range(len(xs)):
            self.word_starts.append(i)
            i += len(xs[n]) + 1

    def find(self, y):
        """Returns starting position of the substring y in the string used for
        building the Suffix tree.

        :param y: String
        :return: Index of the starting position of string y in the string used for building the Suffix tree
                 -1 if y is not a substring.
        """
        node = self.root
        while True:
            edge = self._edgeLabel(node, node.parent)
            if edge.startswith(y):
                return node.idx

            i = 0
            while (i < len(edge) and edge[i] == y[0]):
                y = y[1:]
                i += 1

            if i != 0:
                if i == len(edge) and y != '':
                    pass
                else:
                    return -1

            node = node._get_transition_link(y[0])
            if not node:
                return -1

    def find_all(self, y):
        y_input = y
        node = self.root
        while True:
            edge = self._edgeLabel(node, node.parent)
            if edge.startswith(y):
                break

            i = 0
            while (i < len(edge) and edge[i] == y[0]):
                y = y[1:]
                i += 1

            if i != 0:
                if i == len(edge) and y != '':
                    pass
                else:
                    return []

            node = node._get_transition_link(y[0])
            if not node:
                return []

        leaves = node._get_leaves()
        return [n.idx for n in leaves]

    def _edgeLabel(self, node, parent):
        """Helper method, returns the edge label between a node and it's parent"""
        return self.word[node.idx + parent.depth: node.idx + node.depth]

    def _terminalSymbolsGenerator(self):
        """Generator of unique terminal symbols used for building the Generalized Suffix Tree.
        Unicode Private Use Area U+E000..U+F8FF is used to ensure that terminal symbols
        are not part of the input string.
        """
        if self.end_symbol_set == "unicode":
            symbol_set = UNICODE_SYMBOLS
        else:
            symbol_set = ASCII_SYMBOLS

        for c in symbol_set:
            yield (c)
        raise ValueError("To many input strings.")


class _SNode():
    """Class representing a Node in the Suffix tree."""

    def __init__(self, idx=-1, parentNode=None, depth=-1):
        # Links
        self._suffix_link = None
        self.transition_links = []
        # Properties
        self.idx = idx
        self.depth = depth
        self.parent = parentNode
        self.generalized_idxs = {}

    def __str__(self):
        return ("SNode: idx:" + str(self.idx) + " depth:" + str(self.depth) +
                " transitons:" + str(self.transition_links))

    def _add_suffix_link(self, snode):
        self._suffix_link = snode

    def _get_suffix_link(self):
        if self._suffix_link != None:
            return self._suffix_link
        else:
            return False

    def _get_transition_link(self, suffix):
        for node, _suffix in self.transition_links:
            if _suffix == '__@__' or suffix == _suffix:
                return node
        return False

    def _add_transition_link(self, snode, suffix=''):
        tl = self._get_transition_link(suffix)
        if tl:  # TODO: imporve this.
            self.transition_links.remove((tl, suffix))
        self.transition_links.append((snode, suffix))

    def _has_transition(self, suffix):
        for node, _suffix in self.transition_links:
            if _suffix == '__@__' or suffix == _suffix:
                return True
        return False

    def is_leaf(self):
        return self.transition_links == []

    def _traverse(self, f):
        for (node, _) in self.transition_links:
            node._traverse(f)
        f(self)

    def _get_leaves(self):
        if self.is_leaf():
            return [self]
        else:
            return [x for (n, _) in self.transition_links for x in n._get_leaves()]
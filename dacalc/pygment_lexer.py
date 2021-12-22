from pygments.lexer import RegexLexer

import pygments.token as t

__all__ = ["DALexer"]


class DALexer(RegexLexer):
    '''
    Stripped down lexer of Dimensional Analysis scripts
    (only used for syntax highlighting)
    '''

    name = 'DA'
    aliases = ['da', 'DACalc', 'DACalculator', 'dacalc']
    filenames = ['*.da']

    tokens = {
        'root': [
            (r'\s+', t.Text),
            (r'\#[^\n]*\n?', t.Comment),
            (r'\[[^\]]*\]?', t.String),
            (r'"[^"]*"?', t.Generic.Strong),
            (r'(def|use|output|dim|import|image|analyze|\?)', t.Keyword),
            ((r'(sqrt|sin|cos|tan|asin|acos|atan|atan2|'
              r'abs|log|ln|log2|log10|pow)'), t.Name.Function),
            (r'[a-zA-Z_]+', t.Name.Variable),
            (r'(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?', t.Number),
            (r'[\+\-\*/,=\(\)\{\}\^\%]+', t.Operator)
        ]
    }

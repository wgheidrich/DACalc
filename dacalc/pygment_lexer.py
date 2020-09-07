from pygments.lexer import RegexLexer
from pygments.token import *


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
            (r'\s+', Text),
            (r'\#[^\n]*\n', Comment),
            (r'\[[^\]]*\]', String),
            (r'"[^"]*"', Generic.Strong),
            (r'(def|use|import|analyze|\?)', Keyword),
            (r'(sqrt|sin|cos|tan|asin|acos|atan|atan2|abs|log|ln|log2|log10|pow)', Name.Function),
            (r'[a-zA-Z_]+', Name.Variable),
            (r'(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?', Number),
            (r'[\+\-\*/,=\(\)\^\%]+', Operator)
        ]
    }


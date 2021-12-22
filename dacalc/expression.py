# math expression tree classes
# Author: Wolfgang Heidrich

class Expression:
    '''
    Base class for math expressions -- needs to be subclassed
    '''

    def __init__(self):
        'initialization'
        self.children = []
        pass

    def _eval_children(self, vars):
        'evaluate all child nodes'
        new_children = [c.eval(vars) if isinstance(c, Expression) else c
                        for c in self.children]
        return new_children

    def _is_fully_defined(self, children):
        '''
        check if all child nodes are fully defined
        (i.e. don't have free variables)
        this assumes that eval() has already been run on the tree
        to simplify it as much as possible
        '''
        for c in children:
            if isinstance(c, Expression):
                return False
        return True

    def _format_children(self, priority, readable):
        '''
        list of strings representing the children of this node
        '''
        return [c.format(priority, readable)
                if isinstance(c, Expression) else str(c)
                for c in self.children]

    def get_undefined(self):
        '''
        return dict of undefined variables by recursively traversing the tree

        Returns:
            dict of variable names (with value None)
        '''
        undefined = {}
        for c in self.children:
            if isinstance(c, Expression):
                undefined.update(c.get_undefined())
        return undefined

    def eval(self, vars):
        '''
        Evaluate parse tree

        Must be implemented by subclasses

        Parameters:
            vars :      dictionary of variables and their values

        Returns:
            EITHER a value (int, float,..) if all variables used in the
                   expression are provided
            OR     a new (simplified) parse tree with the provided variables
                   subsituted for their values, and then simplified
        '''
        raise(NotImplementedError)

    def format(self, priority, readable=True):
        '''
        format the expression as a string

        must be implemented by subclasses

        Parameters:
            priority : priority level of the parent operator (determines the
                       use of brackets

            readable : produce readable string vs. internal representation

        Returns:
            string representation
        '''
        raise(NotImplementedError)

    def __str__(self):
        '''
        string conversion
        '''
        return self.format(0)

    def __repr__(self):
        '''
        representation

        same as string but always fully bracketed
        '''
        return self.format(0, False)


class VarExpr(Expression):
    '''
    Expression node for a variable
    '''

    def __init__(self, name):
        '''
        Initialization

        Parameters:
            name :      the name of the variable
        '''
        super().__init__()
        self.name = name

    def eval(self, vars):
        '''
        Evaluate variable

        If the the name of this variable is in the provided variable list,
        it is substituted for its value and returned. Otherwise, self is
        returned unaltered.

        Parameters:
            vars :      dictionary of variables and their values

        Returns:
            EITHER value OR self
        '''
        val = vars.get(self.name, self)
        if isinstance(val, Expression):
            # recursively try to substitute variables
            if not isinstance(val, VarExpr) or val.name != self.name:
                val = val.eval(vars)
        return val

    def get_undefined(self):
        '''
        return dict of undefined variable (i.e just the name of self)

        Returns:
            dictionary with the name of this variable
        '''
        return {self.name: None}

    def format(self, priority, readable=True):
        '''
        format the expression as a string

        Parameters:
            priority : ignored for VarExpr

            readable : produce readable string vs. internal representation
                       (ignored for VarExpr)

        Returns:
            variable name as string
        '''
        return self.name


class UniOpExpr(Expression):
    '''
    Unary Operator Expression

    one of the operators - or \u221a
    (where \u221aa denotes the square root of a)
    '''

    op_priorities = {
        '-': 3,
        '\u221a': 7
    }

    def __init__(self, op, operand):
        '''
        Initialization

        Parameters:
            op :        string of the operator symbol (see list above)

            operand :   operand (node of class Expression or value)
        '''
        super().__init__()
        if op in UniOpExpr.op_priorities.keys():
            self.op = op
        else:
            raise NotImplementedError("Unknown operator " + op)
        self.children.append(operand)

    def eval(self, vars):
        '''
        evaluate unary operator

        recursively evaluate the operand and apply operator

        Parameters:
            vars :     dictionary of variables and their values

        Returns:
            either a value, or new Expression if the operand is not
            fully defined
        '''
        children = self._eval_children(vars)
        if self._is_fully_defined(children):
            operand = children[0]
            if self.op == '-':
                return - operand
            elif self.op == '\u221a':
                # !!! should be cn.sqrt -- fix later
                return operand ** 0.5
            else:
                # we shouldn't get here, since this is checked in __init__
                raise NotImplementedError("Unknown operator " + self.op)
        else:
            return UniOpExpr(self.op, children[0])

    def format(self, priority, readable=True):
        '''
        format the expression as a string

        Parameters:
            priority :  operator priority of the parent
                        the output will be bracketed if this is lower than the
                        priority of self or if readable is False

            readable :  produce readable string vs. internal representation

        Returns:
            operator expression as string
        '''
        my_priority = UniOpExpr.op_priorities[self.op]
        child_strs = self._format_children(my_priority, readable)
        if not readable or priority > my_priority:
            return '(' + self.op + child_strs[0] + ')'
        else:
            return self.op + child_strs[0]


class BinOpExpr(Expression):
    '''
    Binary Operator Expression

    one of the operators +,-,*,/,^, or %
    (where % denotes the root, as in 3%a is the cube root of a)
    '''

    op_priorities = {
        '+': 0,
        '-': 1,  # higher than + so that a-(b+c) gets bracketed the right way
        '*': 3,
        '/': 4,  # same reason as -
        '^': 5,
        '\u221a': 6
    }

    def __init__(self, op, operand1, operand2):
        '''
        Initialization

        Parameters:
            op :        string of the operator name (one of the ones listed

            operand1 :  first operand (node of class Expression or value)

            operand2 :  second operand (node of class Expression or value)
        '''
        super().__init__()
        if op in BinOpExpr.op_priorities.keys():
            self.op = op
        else:
            raise NotImplementedError("Unknown operator " + op)
        self.children.append(operand1)
        self.children.append(operand2)

    def eval(self, vars):
        '''
        evaluate binary operator

        recursively evaluate the operands and apply operator

        Parameters:
            vars :      dictionary of variables and their values

        Returns:
            either a value, or new Expression if the operands are not
            fully defined
        '''
        children = self._eval_children(vars)
        if self._is_fully_defined(children):
            op1 = children[0]
            op2 = children[1]
            if self.op == '+':
                return op1 + op2
            elif self.op == '-':
                return op1 - op2
            elif self.op == '*':
                return op1 * op2
            elif self.op == '/':
                return op1 / op2
            elif self.op == '^':
                return op1 ** op2
            elif self.op == '\u221a':
                # !!! should be cn.root -- fix later
                return op2 ** (1 / op1)
            else:
                # we shouldn't get here, since this is checked in __init__
                raise NotImplementedError("Unknown operator " + self.op)
        else:
            return BinOpExpr(self.op, *children)

    def format(self, priority, readable=True):
        '''
        format the expression as a string

        Parameters:
            priority :  operator priority of the parent
                        the output will be bracketed if this is lower than the
                        priority of self

            readable :  produce readable string vs. internal representation

        Returns:
            operator expression as string
        '''
        my_priority = BinOpExpr.op_priorities[self.op]
        child_strs = self._format_children(my_priority, readable)
        if not readable or priority > my_priority:
            return '(' + child_strs[0] + self.op + child_strs[1] + ')'
        else:
            return child_strs[0] + self.op + child_strs[1]


class FuncExpr(Expression):
    '''
    Function Expression with any number of arguments
    '''
    def __init__(self, name, func, *argv):
        '''
        Initialization

        Parameters:
            name :      name of the function for string representations

            func :      the python function to call

            *argv :     sequence of function arguments
                        (either numerical values or nodes of type Expression)
        '''
        super().__init__()
        self.name = name
        self.func = func
        self.children = list(argv)

    def eval(self, vars):
        '''
        evaluate binary operator

        recursively evaluate the operands and apply operator

        Parameters:
            vars :      dictionary of variables and their values

        Returns:
            either a value, or self if the operand is not fully defined
        '''
        children = self._eval_children(vars)
        if self._is_fully_defined(children):
            return self.func(*children)
        else:
            return FuncExpr(self.name, self.func, *children)

    def format(self, priority, readable=True):
        '''
        format the expression as a string

        Parameters:
            priority :  ignored for BinFuncExpr

            readable :  produce readable string vs. internal representation

        Returns:
            string with function call syntax
        '''
        child_strs = self._format_children(priority, readable)
        if readable:
            return self.name + '(' + ', '.join(child_strs) + ')'
        else:
            return str(self.func) + '(' + ', '.join(child_strs) + ')'


def _test():
    import math as m

    expr = BinOpExpr('+',
                     BinOpExpr('-',
                               1,
                               FuncExpr("cos",
                                        m.cos,
                                        BinOpExpr('/', m.pi, 3))),
                     BinOpExpr('*',
                               VarExpr("x"),
                               VarExpr("y")))

    title = "Expression Tree Example"
    print('\n' + title + '\n' + '=' * len(title) + '\n')
    print("Undefined vars:\t", expr.get_undefined())
    print("Internal repr.:\t", expr.__repr__())
    print("Readable repr:\t", expr)
    eval = expr.eval({})
    print("Simplified:\t", eval)
    evalx = eval.eval({"x": 1})
    print("x = 1:\t\t", evalx)
    evaly = eval.eval({"x": 1, "y": 2})
    print("x = 1, y = 2:\t", evaly)
    print()


if __name__ == '__main__':
    _test()

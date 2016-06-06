#include <stdio.h>
#include "Python.h"

static PyObject *keywdarg_parrot(self, args, keywds)
    PyObject *self;
    PyObject *args;
    PyObject *keywds;
    {  
        int voltage;
        char *state = "a stiff";
        char *action = "voom";
        char *type = "Norwegian Blue";
        static char *kwlist[] = {"voltage", "state", "action", "type", NULL};

        if (!PyArg_ParseTupleAndKeywords(args, keywds, "i|sss", kwlist, &voltage, &state, &action, &type))
            return NULL;

        printf("-- This parrot wouldn't %s if you put %i Volts through it.\n", action, voltage);
        printf("-- Lovely plumage, the %s -- It's %s!\n", type, state);

        Py_INCREF(Py_None);

        return Py_None;
    }

static PyMethodDef keywdarg_methods[] = {
/*  The cast of the function is necessary since PyCFunction values
    only take two PyObject* parameters, and keywdarg_parrot() takes three. */
    {"parrot", (PyCFunction)keywdarg_parrot, METH_VARARGS|METH_KEYWORDS},
    {NULL,  NULL}   /* sentinel */
};

void init_keywdarg(){
/* Create the module and add the functions */
    Py_InitModule("_keywdarg", keywdarg_methods);
}


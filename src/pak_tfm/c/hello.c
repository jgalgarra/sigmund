#include <Python.h>
#include <math.h>

struct module_state {
    PyObject *error;
};

/*
static PyObject* helloworld(PyObject* self, PyObject* args)
{
	PyObject * myobj = NULL;
	// printf("%f\n",&args);
	if (!PyArg_parsetuple(args,"f", &myobj))
        return NULL;
    printf("%f\n",myobj);
    return Py_BuildValue("f", myobj);
}
*/

static PyObject* helloworld(PyObject* self, PyObject * args)
{
	//return Py_BuildValue("s", "Hello, Python extensions!!");
	//float* myobj,voltage;
	// if (!PyArg_ParseTuple(args, "f", &myobj)) {
    //    return(NULL);}
    PyObject *num = PyNumber_Float(PyTuple_GetItem(args, 0));
	//voltage = PyFloat_AsDouble(num);
	//printf("voltage %f ", voltage);
    return Py_BuildValue("f", expm1f(PyFloat_AsDouble(num)));
}

static char helloworld_docs[] =
    "helloworld( ): Any message you want to put here!!\n";

static PyMethodDef helloworld_methods[] = {
    {"helloworld", (PyCFunction)helloworld, 
     METH_VARARGS, helloworld_docs},
    {NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "helloworld",
        NULL,
        sizeof(struct module_state),
        helloworld_methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject * PyInit_helloworld(void)
/*
void inithelloworld(void)*/
{
    PyObject *module = PyModule_Create(&moduledef);
}
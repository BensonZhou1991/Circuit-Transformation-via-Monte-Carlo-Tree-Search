OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
cx q[2],q[1];
cx q[1],q[2];
cx q[3],q[2];
cx q[2],q[3];
cx q[4],q[3];
cx q[3],q[4];
cx q[4],q[1];
cx q[0],q[4];
cx q[1],q[0];
cx q[1],q[4];
cx q[0],q[4];
cx q[1],q[0];
cx q[4],q[1];
cx q[0],q[4];
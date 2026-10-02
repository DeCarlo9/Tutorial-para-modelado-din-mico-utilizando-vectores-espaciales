import time
import matplotlib.pyplot as plt
import numpy as np
import pybullet as p
import pybullet_data


def cinematica_directa(q):
    q1, q2, q3 = q[0], q[1], q[2]

    l1 = 0.4
    l2 = 1.0
    l3 = 1.0

    x = np.cos(q1) * (l3 * np.cos(q2 + q3) + l2 * np.cos(q2))
    y = np.sin(q1) * (l3 * np.cos(q2 + q3) + l2 * np.cos(q2))
    z = l1 + l3 * np.sin(q2 + q3) + l2 * np.sin(q2)

    pos = np.array([x, y, z])

    return pos


def cinematica_inversa(t):
    """Calcula qd y qdp para describir una trayectoria circular en el plano X-Y."""
    l1 = 0.4
    l2 = 1.0
    l3 = 1.0

    # Parámetros de la trayectoria circular
    xc = 1.2
    yc = 0.8
    zc = 0.5
    R = 0.25
    f = 0.2
    omega = 2 * np.pi *f

    # Posición y velocidad deseada en el espacio cartesiano
    x_d = np.array([xc + R * np.cos(omega * t), yc + R * np.sin(omega * t), zc])
    xdp_d = np.array([-R * omega * np.sin(omega * t), R * omega * np.cos(omega * t), 0.0])

    x, y, z = x_d[0], x_d[1], x_d[2]

    # Cinemática inversa de posición
    q1 = np.arctan2(y, x)
    r = np.sqrt(x**2 + y**2)
    zp = z - l1

    D = (r**2 + zp**2 - l2**2 - l3**2) / (2 * l2 * l3)
    #D = np.clip(D, -1.0, 1.0)

    s3 = -np.sqrt(1 - D**2)  # configuración codo arriba
    q3 = np.arctan2(s3, D)
    q2 = np.arctan2(zp, r) - np.arctan2(l3 * s3, l2 + l3 * D)

    qd = np.array([q1, q2, q3])

    # Jacobiano analítico en qd
    rp_q2p = -l2 * np.sin(q2) - l3 * np.sin(q2 + q3)
    rp_q3p = -l3 * np.sin(q2 + q3)

    J11 = -r * np.sin(q1)
    J12 = rp_q2p * np.cos(q1)
    J13 = rp_q3p * np.cos(q1)

    J21 = r * np.cos(q1)
    J22 = rp_q2p * np.sin(q1)
    J23 = rp_q3p * np.sin(q1)

    J31 = 0.0
    J32 = l2 * np.cos(q2) + l3 * np.cos(q2 + q3)
    J33 = l3 * np.cos(q2 + q3)

    J = np.array([[J11, J12, J13], [J21, J22, J23], [J31, J32, J33]])

    # Cinemática inversa de velocidad
    qdp = np.linalg.solve(J, xdp_d)

    return qd, qdp


def modelo_antropomorfico(q, qp, tau, F_ext):
    # Parámetros físicos
    l2, lc2 = 1.0, 0.5
    l3, lc3 = 1.0, 0.5
    m1, m2, m3 = 1.0, 1.0, 1.0
    I1, I2, I3 = 0.001, 0.001, 0.001
    g = 9.81

    q1, q2, q3 = q[0], q[1], q[2]
    q1p, q2p, q3p = qp[0], qp[1], qp[2]

    # Matriz de inercias
    M11 = (
        I1
        + I2
        + I3
        + 0.5 * l2**2 * m3
        + 0.5 * lc2**2 * m2
        + 0.5 * lc3**2 * m3
        + 0.5 * l2**2 * m3 * np.cos(2 * q2)
        + 0.5 * lc2**2 * m2 * np.cos(2 * q2)
        + 0.5 * lc3**2 * m3 * np.cos(2 * q2 + 2 * q3)
        + l2 * lc3 * m3 * np.cos(q3)
        + l2 * lc3 * m3 * np.cos(2 * q2 + q3)
    )
    M22 = (
        m3 * l2**2
        + 2 * m3 * np.cos(q3) * l2 * lc3
        + m2 * lc2**2
        + m3 * lc3**2
        + I2
        + I3
    )
    M23 = m3 * lc3**2 + l2 * m3 * np.cos(q3) * lc3 + I3
    M33 = m3 * lc3**2 + I3

    M = np.array([[M11, 0.0, 0.0], [0.0, M22, M23], [0.0, M23, M33]])

    # Matriz de coriolis
    term1 = (
        0.5 * m3 * np.sin(2 * q2) * l2**2
        + m3 * np.sin(2 * q2 + q3) * l2 * lc3
        + 0.5 * m2 * np.sin(2 * q2) * lc2**2
        + 0.5 * m3 * np.sin(2 * q2 + 2 * q3) * lc3**2
    )
    term2 = 0.5 * lc3 * m3 * (
        lc3 * np.sin(2 * q2 + 2 * q3)
        + l2 * np.sin(q3)
        + l2 * np.sin(2 * q2 + q3)
    )

    C11 = -q2p * term1 - q3p * term2
    C12 = -q1p * term1
    C13 = -q1p * term2
    C21 = q1p * term1
    C22 = -l2 * lc3 * m3 * q3p * np.sin(q3)
    C23 = -l2 * lc3 * m3 * np.sin(q3) * (q2p + q3p)
    C31 = q1p * term2
    C32 = l2 * lc3 * m3 * q2p * np.sin(q3)
    C33 = 0.0

    C = np.array([[C11, C12, C13], [C21, C22, C23], [C31, C32, C33]])

    # Vector de gravedad
    g1 = 0.0
    g2 = g * m3 * (lc3 * np.cos(q2 + q3) + l2 * np.cos(q2)) + g * lc2 * m2 * np.cos(q2)
    g3 = g * lc3 * m3 * np.cos(q2 + q3)
    G = np.array([g1, g2, g3])

    # Jacobiano para aplicar la perturbación externa 
    r = l2 * np.cos(q2) + l3 * np.cos(q2 + q3)
    rp_q2p = -l2 * np.sin(q2) - l3 * np.sin(q2 + q3)
    rp_q3p = -l3 * np.sin(q2 + q3)

    J = np.array([
        [-r * np.sin(q1), rp_q2p * np.cos(q1), rp_q3p * np.cos(q1)],
        [r * np.cos(q1), rp_q2p * np.sin(q1), rp_q3p * np.sin(q1)],
        [0.0, l2 * np.cos(q2) + l3 * np.cos(q2 + q3), l3 * np.cos(q2 + q3)],
    ])

    tau_ext = J.T @ F_ext
    qpp = np.linalg.solve(M, tau + tau_ext - C @ qp - G)
    return qpp


# --- Inicialización de PyBullet ---
physicsClient = p.connect(p.GUI)
p.setAdditionalSearchPath(pybullet_data.getDataPath())
p.setGravity(0, 0, -9.81)

robot_id = p.loadURDF("robot_3gdl.urdf", [0, 0, 0], useFixedBase=True)
joint_indices = [0, 1, 2]

for joint in joint_indices:
    p.setJointMotorControl2(
        bodyUniqueId=robot_id,
        jointIndex=joint,
        controlMode=p.VELOCITY_CONTROL,
        force=0,
    )
    p.changeDynamics(
        bodyUniqueId=robot_id,
        linkIndex=joint,
        jointDamping=0.0,
        linearDamping=0.0,
        angularDamping=0.0,
        lateralFriction=0.0,
        spinningFriction=0.0,
        rollingFriction=0.0,
    )

p.resetDebugVisualizerCamera(
    cameraDistance=3.0,
    cameraYaw=90,
    cameraPitch=-60,
    cameraTargetPosition=[1, 0, 1],
)

# Fuerza externa en coordenadas del mundo
F_ext = np.array([0.0, 0.0, -50.0])

#Ganancias
Kp = np.diag([100.0, 100.0, 100.0])
Kd = np.diag([15.0, 15.0, 15.0])
Ki = np.diag([0.0, 0.0, 0.0])

t_experimento = 10.0
dt = 0.001
p.setTimeStep(dt)
pasos = int(t_experimento / dt)

# Asignación de condiciones iniciales
q_init = np.array([np.pi/6, np.pi/4, -np.pi/3])

# Restablecer robot en PyBullet
for i, j_idx in enumerate(joint_indices):
    p.resetJointState(robot_id, j_idx, q_init[i])

# Asignar estado al modelo analítico
q_analitico = q_init.copy()
qp_analitico = np.array([0.0, 0.0, 0.0])

# Acumuladores del error integral
e_int_pb = np.array([0.0, 0.0, 0.0])
e_int_sim = np.array([0.0, 0.0, 0.0])

contenedor_t = []
contenedor_q_bullet = []
contenedor_q_analitico = []
contenedor_qd = []
contenedor_pos_analitica = []
contenedor_pos_pybullet = []

for paso in range(pasos):
    t = paso * dt

    # Obtener trayectoria deseada por cinemática inversa
    qd, qdp = cinematica_inversa(t)

    # Estados PyBullet
    q_pb = np.array([p.getJointState(robot_id, j)[0] for j in joint_indices])
    qp_pb = np.array([p.getJointState(robot_id, j)[1] for j in joint_indices])

    # Obtener la posición de la punta en coordenadas del mundo
    link_state_end = p.getLinkState(robot_id, 2)
    pos_origen_link2 = np.array(link_state_end[0])
    R_matrix_end = np.array(p.getMatrixFromQuaternion(link_state_end[1])).reshape(3, 3)
    pos_pybullet = pos_origen_link2 + (R_matrix_end @ np.array([0.5, 0.0, 0.0]))

    # Aplicar la fuerza externa en la punta
    p.applyExternalForce(
        objectUniqueId=robot_id,
        linkIndex=2,
        forceObj=F_ext,          # Fuerza en el mundo
        posObj=pos_pybullet,     # Punto de aplicación exacto en el mundo
        flags=p.WORLD_FRAME,
    )

    # --- Control PID (PyBullet) ---
    e_pb = qd - q_pb
    ep_pb = qdp - qp_pb
    e_int_pb += e_pb * dt

    tau_pb = Kp @ e_pb + Kd @ ep_pb + Ki @ e_int_pb

    p.setJointMotorControlArray(
        bodyUniqueId=robot_id,
        jointIndices=joint_indices,
        controlMode=p.TORQUE_CONTROL,
        forces=tau_pb,
    )
    p.stepSimulation()

    # --- Control PID (Modelo Analítico) ---
    e_sim = qd - q_analitico
    ep_sim = qdp - qp_analitico
    e_int_sim += e_sim * dt

    tau_sim = Kp @ e_sim + Kd @ ep_sim + Ki @ e_int_sim

    qpp_sim = modelo_antropomorfico(q_analitico, qp_analitico, tau_sim, F_ext)
    qp_analitico += qpp_sim * dt
    q_analitico += qp_analitico * dt

    pos_analitica = cinematica_directa(q_analitico)

    # Guardar variables
    contenedor_t.append(t)
    contenedor_q_bullet.append(q_pb)
    contenedor_q_analitico.append(q_analitico.copy())
    contenedor_qd.append(qd)
    contenedor_pos_analitica.append(pos_analitica)
    contenedor_pos_pybullet.append(pos_pybullet)

p.disconnect()

# --- Gráficas de Resultados ---
q_bullet_deg = np.array(contenedor_q_bullet) * (180.0 / np.pi)
q_analitico_deg = np.array(contenedor_q_analitico) * (180.0 / np.pi)
qd_deg = np.array(contenedor_qd) * (180.0 / np.pi)

fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
nombres = ["q1", "q2", "q3"]

for i in range(3):
    axs[i].plot(contenedor_t, q_bullet_deg[:, i], "b-", label=f"{nombres[i]} PyBullet", linewidth=2)
    axs[i].plot(contenedor_t, q_analitico_deg[:, i], "r--", label=f"{nombres[i]} Analítico", linewidth=2)
    axs[i].plot(contenedor_t, qd_deg[:, i], "k:", label="qd (Cin. Inversa)")
    axs[i].set_ylabel(f"{nombres[i]} [°]")
    axs[i].grid(True)
    axs[i].legend(loc="upper right")

axs[2].set_xlabel("Tiempo [s]")
fig.suptitle("Analítico vs Pybullet", fontsize=14)
plt.tight_layout()
plt.show()

# --- Exportar Datos a Archivo de Texto ---
t_arr = np.array(contenedor_t).reshape(-1, 1)
q_bullet_deg = np.array(contenedor_q_bullet)
q_analitico_deg = np.array(contenedor_q_analitico)
qd_deg = np.array(contenedor_qd)
pos_analitico = np.array(contenedor_pos_analitica)
pos_pb = np.array(contenedor_pos_pybullet)

datos_completos = np.hstack((t_arr, q_bullet_deg, q_analitico_deg, qd_deg, pos_analitico, pos_pb))

nombre_archivo = "datos.txt"
np.savetxt(
    nombre_archivo,
    datos_completos,
    delimiter=",",
    fmt="%.6f",
    comments=""
)

print(f"Datos guardados en '{nombre_archivo}'")

% Calculate time derivative of joint inertia matrix for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m2DE_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'palh2m2DE_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:29
% DurationCPUTime: 0.36s
% Computational Cost: add. (117->68), mult. (294->123), div. (0->0), fcn. (78->6), ass. (0->41)
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t4 = t12 * rSges(7,1) + t15 * rSges(7,2);
t59 = qJD(4) * t4;
t13 = sin(qJ(3));
t26 = rSges(7,1) * t15 - rSges(7,2) * t12;
t2 = t26 * qJD(4);
t16 = cos(qJ(3));
t36 = qJD(3) * t16;
t58 = t13 * t2 + t4 * t36;
t14 = sin(qJ(2));
t17 = cos(qJ(2));
t38 = qJD(2) * t17;
t57 = (t14 * t2 + t4 * t38) * m(7);
t50 = (m(5) * rSges(5,2));
t55 = -2 * rSges(5,1) * t50 + 2 * Icges(5,4);
t54 = -2 * m(3) * rSges(3,1) * rSges(3,2) + 2 * Icges(3,4);
t53 = pkin(5) * m(7);
t52 = (m(5) + m(6));
t51 = m(6) + m(7);
t49 = rSges(3,3) * m(3);
t48 = rSges(5,3) * m(5);
t47 = rSges(6,3) * m(6);
t18 = m(5) * rSges(5,1);
t42 = pkin(5) * m(6) + t18;
t41 = m(3) * pkin(1);
t34 = pkin(1) + pkin(2);
t10 = pkin(3) + t34;
t40 = qJD(4) * (t17 * pkin(4) + t16 * pkin(5) + t10);
t39 = qJD(2) * t14;
t37 = qJD(3) * t13;
t35 = 0.2e1 * m(7);
t33 = t13 * t50;
t32 = t16 * t50;
t23 = t35 * t59;
t22 = 0.2e1 * m(6) * rSges(6,1) + (t10 + t26) * t35;
t11 = (m(4) + t52);
t6 = t51 * pkin(5) + t18;
t5 = t42 + t53;
t3 = -pkin(4) * t39 - pkin(5) * t37;
t1 = [(t38 * t54 + (0.2e1 * (-t13 * t5 - t32) * qJD(3) - t23) * pkin(4)) * t17 - ((2 * rSges(3,1) * t41) + t14 * t54 + ((2 * t11 * pkin(1)) + (2 * t52 * pkin(2)) + (2 * m(4) * rSges(4,1)) + 0.2e1 * t5 * t16 + t22 - 0.2e1 * t33) * pkin(4)) * t39 + (-pkin(5) * t23 + t36 * t55) * t16 - (t22 * pkin(5) + t13 * t55 + 0.2e1 * t34 * t42) * t37 - 0.2e1 * ((rSges(3,2) * t41) + (((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3)) - Icges(3,1) + Icges(3,2) + (m(7) + t11) * pkin(4) ^ 2) * t14) * t38 - 0.2e1 * ((t34 * t50) + (((rSges(5,1) ^ 2 - rSges(5,2) ^ 2) * m(5)) - Icges(5,1) + Icges(5,2) + t51 * pkin(5) ^ 2) * t13) * t36 - 0.2e1 * m(7) * t10 * t59; ((-rSges(3,1) * t49 + Icges(3,5)) * t17 + (rSges(3,2) * t49 - Icges(3,6)) * t14) * qJD(2) + ((-(rSges(4,3) * m(4)) - t47 - t48) * t38 + t57) * pkin(4); 0; ((-rSges(5,1) * t48 + Icges(5,5)) * t16 - (-rSges(5,2) * t48 + Icges(5,6)) * t13) * qJD(3) + (m(7) * t58 - t36 * t47) * pkin(5); (qJD(2) - qJD(3)) * ((t6 * t13 + t32) * t17 - (t6 * t16 - t33) * t14) * pkin(4); 0; ((-t12 * t3 - t15 * t40) * rSges(7,2) + (-t12 * t40 + t15 * t3) * rSges(7,1)) * m(7); pkin(4) * t57; t58 * t53; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

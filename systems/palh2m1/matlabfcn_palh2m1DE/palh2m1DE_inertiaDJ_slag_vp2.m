% Calculate time derivative of joint inertia matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = palh2m1DE_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1DE_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1DE_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:08
% DurationCPUTime: 0.36s
% Computational Cost: add. (255->77), mult. (302->117), div. (0->0), fcn. (95->22), ass. (0->55)
t34 = qJ(2) + qJ(3);
t20 = sin(t34);
t72 = -0.2e1 * t20;
t32 = -qJ(4) + qJ(2);
t25 = qJ(3) + t32;
t10 = sin(t25);
t27 = -qJD(4) + qJD(2);
t17 = qJD(3) + t27;
t62 = t17 * t10;
t26 = qJD(4) + qJD(2);
t16 = qJD(3) + t26;
t31 = qJ(4) + qJ(2);
t24 = qJ(3) + t31;
t9 = sin(t24);
t64 = t16 * t9;
t71 = (-t62 / 0.2e1 + t64 / 0.2e1) * mrSges(6,1);
t12 = cos(t25);
t61 = t17 * t12;
t11 = cos(t24);
t63 = t16 * t11;
t70 = (t61 + t63) * mrSges(6,2);
t56 = t27 * cos(t32);
t57 = t27 * sin(t32);
t58 = t26 * cos(t31);
t59 = t26 * sin(t31);
t68 = (t56 / 0.2e1 + t58 / 0.2e1) * mrSges(6,2) + (-t57 / 0.2e1 + t59 / 0.2e1) * mrSges(6,1);
t67 = -2 * pkin(1);
t66 = pkin(3) / 0.2e1;
t42 = m(5) + m(6);
t41 = pkin(1) + pkin(4);
t23 = cos(t34);
t29 = qJD(2) + qJD(3);
t60 = t23 * t29;
t39 = cos(qJ(3));
t13 = t39 * pkin(3) + pkin(2);
t36 = sin(qJ(3));
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t55 = qJD(4) * (-t37 * t36 * pkin(3) + t13 * t40 + t41);
t54 = qJD(2) * t40;
t53 = -0.2e1 * qJD(4) * t41;
t48 = m(6) * pkin(4) + mrSges(5,1);
t14 = t42 * pkin(3) + mrSges(4,1);
t45 = Ifges(4,6) * t29 * t20 + (mrSges(5,3) * pkin(3) - Ifges(4,5)) * t60 + pkin(3) * t71 + t66 * t70;
t44 = pkin(3) ^ 2;
t43 = 0.2e1 * qJ(2);
t38 = cos(qJ(4));
t35 = sin(qJ(4));
t33 = t43 + qJ(3);
t30 = m(4) + t42;
t28 = 0.2e1 * qJD(2) + qJD(3);
t15 = 0.2e1 * t34;
t2 = (-mrSges(4,2) * t39 - t14 * t36) * qJD(3) * pkin(2);
t1 = -t13 * qJD(2) * t37 + (-t36 * t54 + (-t36 * t40 - t37 * t39) * qJD(3)) * pkin(3);
t3 = [((pkin(1) * t14 + pkin(3) * t48) * t72 + mrSges(4,2) * t23 * t67 + 0.2e1 * Ifges(4,4) * cos(t15) - (t42 * t44 - Ifges(4,1) + Ifges(4,2)) * sin(t15)) * t29 + ((-sin(t33) * t28 - qJD(3) * t36) * t14 + (-cos(t33) * t28 - qJD(3) * t39) * mrSges(4,2)) * pkin(2) + (0.2e1 * Ifges(3,4) * cos(t43) - 0.2e1 * ((pkin(1) * t30 + t48) * pkin(2) + (mrSges(3,1) * pkin(1))) * t37 - (t30 * pkin(2) ^ 2 - Ifges(3,1) + Ifges(3,2)) * sin(t43) + mrSges(3,2) * t40 * t67) * qJD(2) + (t38 * t53 + (t61 - t63) * pkin(3) + (t56 - t58) * pkin(2)) * mrSges(6,2) + (t35 * t53 + (-t62 - t64) * pkin(3) + (-t57 - t59) * pkin(2)) * mrSges(6,1); (-Ifges(3,5) * t40 + Ifges(3,6) * t37) * qJD(2) + ((mrSges(4,3) + mrSges(5,3)) * t54 + t68) * pkin(2) + t45; 0.2e1 * t2; t45; t2; (t60 * t72 + (t10 * t11 + t12 * t9) * (t16 / 0.2e1 + t17 / 0.2e1)) * t44 * m(6); (-t1 * t35 - t38 * t55) * mrSges(6,2) + (t1 * t38 - t35 * t55) * mrSges(6,1); (-(-t61 / 0.2e1 - t63 / 0.2e1) * mrSges(6,2) + t71) * pkin(3) + t68 * pkin(2); (t70 + (-t62 + t64) * mrSges(6,1)) * t66; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;

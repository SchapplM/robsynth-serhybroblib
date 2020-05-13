% Calculate Gravitation load on the joints for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m1DE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:52:06
% EndTime: 2020-05-02 23:52:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (124->35), mult. (198->55), div. (0->0), fcn. (56->8), ass. (0->23)
t37 = m(5) + m(6);
t40 = m(4) + t37;
t11 = rSges(6,1) * g(1) + rSges(6,2) * g(2);
t12 = rSges(6,1) * g(2) - rSges(6,2) * g(1);
t14 = sin(qJ(4));
t17 = sin(qJ(1));
t18 = cos(qJ(4));
t21 = cos(qJ(1));
t39 = ((t11 * t14 - t12 * t18) * t21 + t17 * (t11 * t18 + t12 * t14)) * m(6);
t15 = sin(qJ(3));
t19 = cos(qJ(3));
t34 = m(4) * rSges(4,1) + pkin(3) * t37;
t36 = m(4) * rSges(4,2);
t6 = rSges(3,2) * m(3) + t15 * t34 + t19 * t36;
t3 = rSges(3,1) * m(3) + t40 * pkin(2) - t15 * t36 + t19 * t34;
t33 = g(1) * t21 + g(2) * t17;
t16 = sin(qJ(2));
t20 = cos(qJ(2));
t30 = t3 * t20 - t6 * t16 + m(2) * rSges(2,1) + rSges(5,1) * m(5) + pkin(4) * m(6) + (m(3) + t40) * pkin(1);
t8 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3);
t4 = g(3) * t34 + t33 * t36;
t2 = -g(3) * t36 + t33 * t34;
t1 = [(t17 * t8 - t30 * t21) * g(2) + (t17 * t30 + t8 * t21) * g(1) + t39, (-g(3) * t6 + t33 * t3) * t16 + (g(3) * t3 + t33 * t6) * t20, (t15 * t2 + t19 * t4) * t20 + (-t15 * t4 + t19 * t2) * t16, t39];
taug = t1(:);

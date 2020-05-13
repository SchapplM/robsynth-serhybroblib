% Calculate Gravitation load on the joints for
% palh2m2DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% m [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = palh2m2DE_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:28
% EndTime: 2020-05-03 01:06:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (97->38), mult. (137->52), div. (0->0), fcn. (32->8), ass. (0->21)
t12 = sin(qJ(1));
t13 = cos(qJ(4));
t16 = cos(qJ(1));
t7 = rSges(7,1) * g(1) + rSges(7,2) * g(2);
t8 = rSges(7,1) * g(2) - rSges(7,2) * g(1);
t9 = sin(qJ(4));
t30 = (t12 * (t13 * t7 + t8 * t9) - (t13 * t8 - t7 * t9) * t16) * m(7);
t28 = m(6) + m(7);
t27 = m(3) * rSges(3,2);
t26 = m(5) * rSges(5,2);
t25 = pkin(1) + pkin(2);
t24 = t28 * pkin(5) + m(5) * rSges(5,1);
t21 = g(1) * t16 + g(2) * t12;
t20 = (m(4) + m(5) + t28) * pkin(4) + m(3) * rSges(3,1);
t10 = sin(qJ(3));
t11 = sin(qJ(2));
t14 = cos(qJ(3));
t15 = cos(qJ(2));
t19 = t10 * t26 + t11 * t27 - t14 * t24 - t15 * t20 - (pkin(3) + t25) * m(7) - t25 * m(5) - (rSges(6,1) + t25) * m(6) - (m(3) + m(4)) * pkin(1) - m(2) * rSges(2,1) - m(4) * rSges(4,1);
t2 = m(2) * rSges(2,2) - m(3) * rSges(3,3) - m(4) * rSges(4,3) - m(5) * rSges(5,3) - m(6) * rSges(6,3);
t1 = [t30 + (t12 * t2 + t19 * t16) * g(2) + (-t12 * t19 + t2 * t16) * g(1), (g(3) * t27 + t20 * t21) * t11 - t15 * (g(3) * t20 - t21 * t27), (g(3) * t26 + t21 * t24) * t10 - (g(3) * t24 - t21 * t26) * t14, t30];
taug = t1(:);

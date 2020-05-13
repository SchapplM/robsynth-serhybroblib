% Calculate Gravitation load on the joints for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbarprisTE_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisTE_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:21
% EndTime: 2020-05-07 09:01:23
% DurationCPUTime: 0.25s
% Computational Cost: add. (104->42), mult. (160->56), div. (3->4), fcn. (4->2), ass. (0->26)
t17 = (qJ(1) ^ 2);
t19 = (pkin(3) ^ 2);
t47 = 2 * qJ(1) * pkin(3) + t17 + t19;
t16 = qJ(1) * t17;
t18 = pkin(3) * t19;
t34 = qJ(1) * t19;
t46 = -4 * t34 - 2 * t16 - t18;
t21 = (pkin(1) ^ 2);
t23 = (pkin(2) ^ 2);
t24 = (pkin(2) * t23);
t42 = -t24 + (t21 - t47) * pkin(2);
t13 = qJ(1) + pkin(3);
t28 = t13 ^ 2;
t12 = t13 * t28;
t41 = -2 * t12;
t39 = -5 * t17;
t38 = (m(3) * g(2));
t14 = (mrSges(2,2) - mrSges(3,1));
t15 = (mrSges(2,1) + mrSges(3,3));
t1 = (g(1) * t14 - t15 * g(2));
t37 = t21 * t1;
t33 = 2 * t12 * (-mrSges(4,1) * g(2) + mrSges(4,2) * g(1));
t32 = -pkin(2) - t13;
t31 = -pkin(2) + t13;
t22 = sqrt(-((pkin(1) + t32) * (pkin(1) + t31) * (pkin(1) - t31) * (pkin(1) - t32)));
t2 = [-(((mrSges(4,1) * t41 + t42 * t15 + (pkin(3) * t24 + ((t39 - t21) * pkin(3) + t46) * pkin(2)) * m(3)) * g(1) + (mrSges(4,2) * t41 + t42 * t14) * g(2)) * t22 + (2 * (-t37 + (-t18 - 3 * t34 - (3 * t17 + t21) * pkin(3) - t16) * t38) * t24) + (((t37 + t47 * t1 + ((t39 + t21) * pkin(3) + t46) * t38) * pkin(2) + t33) * (pkin(1) - t13) * (pkin(1) + t13)) + (((pkin(3) * t38 + t1) * t24 + t33) * t23)) / t22 / pkin(1) / t28 / pkin(2) / 0.2e1];
taug = t2(:);

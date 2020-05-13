% Calculate potential energy for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1DE_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:57
% EndTime: 2020-05-02 23:51:58
% DurationCPUTime: 0.23s
% Computational Cost: add. (98->43), mult. (120->51), div. (0->0), fcn. (24->8), ass. (0->19)
t20 = m(5) + m(6);
t19 = m(4) * rSges(4,2);
t18 = m(4) + t20;
t17 = m(3) + m(2) + t18;
t12 = cos(qJ(3));
t4 = m(4) * rSges(4,1) + pkin(3) * t20;
t9 = sin(qJ(3));
t1 = rSges(3,1) * m(3) + pkin(2) * t18 + t4 * t12 - t19 * t9;
t10 = sin(qJ(2));
t13 = cos(qJ(2));
t2 = rSges(3,2) * m(3) + t12 * t19 + t4 * t9;
t16 = -m(2) * rSges(2,1) - (pkin(1) + rSges(5,1)) * m(5) - (pkin(1) + pkin(4)) * m(6) - t1 * t13 + t2 * t10 + (-m(3) - m(4)) * pkin(1);
t11 = cos(qJ(4));
t8 = sin(qJ(4));
t7 = m(1) + t17;
t6 = -rSges(6,1) * g(2) + rSges(6,2) * g(1);
t5 = rSges(6,1) * g(1) + rSges(6,2) * g(2);
t3 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3);
t14 = (-g(2) * t3 + (-t5 * t11 + t6 * t8) * m(6) + t16 * g(1)) * cos(qJ(1)) + (t3 * g(1) + (t6 * t11 + t5 * t8) * m(6) + t16 * g(2)) * sin(qJ(1)) + (-g(1) * r_base(1) - g(2) * r_base(2)) * t7 + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (t2 * t13 + t1 * t10 - t7 * r_base(3) - (pkin(6) + rSges(6,3)) * m(6) + rSges(5,2) * m(5) - m(2) * rSges(2,3) - m(1) * rSges(1,3) - t17 * pkin(5)) * g(3);
U = t14;

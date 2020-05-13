% Calculate potential energy for
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m1DE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1DE_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'palh2m1DE_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-02 23:51:58
% EndTime: 2020-05-02 23:51:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (80->40), mult. (114->47), div. (0->0), fcn. (24->8), ass. (0->17)
t38 = m(5) + m(6);
t37 = m(4) * rSges(4,2);
t36 = m(4) + t38;
t24 = m(4) * rSges(4,1) + t38 * pkin(3);
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t21 = rSges(3,1) * m(3) + t36 * pkin(2) + t24 * t31 - t28 * t37;
t22 = rSges(3,2) * m(3) + t24 * t28 + t31 * t37;
t29 = sin(qJ(2));
t32 = cos(qJ(2));
t35 = -m(2) * rSges(2,1) - (pkin(1) + rSges(5,1)) * m(5) - (pkin(1) + pkin(4)) * m(6) - t21 * t32 + t22 * t29 + (-m(3) - m(4)) * pkin(1);
t30 = cos(qJ(4));
t27 = sin(qJ(4));
t26 = -rSges(6,1) * g(2) + rSges(6,2) * g(1);
t25 = rSges(6,1) * g(1) + rSges(6,2) * g(2);
t23 = m(2) * rSges(2,2) + m(3) * rSges(3,3) + m(4) * rSges(4,3) + m(5) * rSges(5,3);
t1 = (-g(2) * t23 + (-t25 * t30 + t26 * t27) * m(6) + t35 * g(1)) * cos(qJ(1)) + (t23 * g(1) + (t25 * t27 + t26 * t30) * m(6) + t35 * g(2)) * sin(qJ(1)) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (t22 * t32 + t21 * t29 - (pkin(6) + rSges(6,3)) * m(6) + rSges(5,2) * m(5) - m(2) * rSges(2,3) - m(1) * rSges(1,3) - (m(3) + m(2) + t36) * pkin(5)) * g(3);
U = t1;

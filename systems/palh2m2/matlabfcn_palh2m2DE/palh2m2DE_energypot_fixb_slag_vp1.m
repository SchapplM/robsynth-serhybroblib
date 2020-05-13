% Calculate potential energy for
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
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = palh2m2DE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'palh2m2DE_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'palh2m2DE_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 01:06:21
% EndTime: 2020-05-03 01:06:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (78->42), mult. (108->50), div. (0->0), fcn. (18->8), ass. (0->15)
t36 = m(6) + m(7);
t27 = pkin(1) + pkin(2);
t21 = m(3) * rSges(3,1) + (m(4) + m(5) + t36) * pkin(4);
t22 = t36 * pkin(5) + m(5) * rSges(5,1);
t30 = sin(qJ(3));
t31 = sin(qJ(2));
t33 = cos(qJ(3));
t34 = cos(qJ(2));
t35 = -m(2) * rSges(2,1) - (rSges(4,1) + pkin(1)) * m(4) - (rSges(6,1) + t27) * m(6) - (pkin(3) + t27) * m(7) - t21 * t34 - t22 * t33 + (rSges(5,2) * t30 - t27) * m(5) + (rSges(3,2) * t31 - pkin(1)) * m(3);
t32 = cos(qJ(4));
t29 = sin(qJ(4));
t24 = -rSges(7,1) * g(2) + rSges(7,2) * g(1);
t23 = rSges(7,1) * g(1) + rSges(7,2) * g(2);
t20 = m(2) * rSges(2,2) - m(3) * rSges(3,3) - m(4) * rSges(4,3) - m(5) * rSges(5,3) - m(6) * rSges(6,3);
t1 = (-g(2) * t20 + (-t23 * t32 + t24 * t29) * m(7) + t35 * g(1)) * cos(qJ(1)) + (t20 * g(1) + (t23 * t29 + t24 * t32) * m(7) + t35 * g(2)) * sin(qJ(1)) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2)) * m(1) + (-m(3) * t34 * rSges(3,2) - m(5) * t33 * rSges(5,2) - m(1) * rSges(1,3) - m(2) * rSges(2,3) - m(4) * rSges(4,2) - m(6) * rSges(6,2) - m(7) * rSges(7,3) - t21 * t31 - t22 * t30) * g(3);
U = t1;

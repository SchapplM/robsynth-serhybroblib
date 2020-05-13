% Calculate potential energy for
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbarprisTE_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energypot_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisTE_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energypot_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisTE_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:06
% EndTime: 2020-05-07 09:01:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (80->43), mult. (124->59), div. (4->3), fcn. (4->2), ass. (0->14)
t36 = (rSges(2,1) * m(2) + rSges(3,3) * m(3));
t53 = (m(3) * qJ(1) + t36);
t52 = 2 * pkin(3);
t35 = (-rSges(2,2) * m(2) + rSges(3,1) * m(3));
t51 = (g(2) * t35);
t49 = (pkin(2) - pkin(3)) * (pkin(2) + pkin(3));
t48 = qJ(1) + pkin(3);
t47 = 2 * ((-rSges(2,3) * m(2) + rSges(3,2) * m(3)) * g(3) + (-rSges(1,1) * g(1) - rSges(1,2) * g(2) - rSges(1,3) * g(3)) * m(1)) * pkin(1);
t46 = -pkin(2) - t48;
t45 = -pkin(2) + t48;
t42 = pkin(1) ^ 2;
t41 = qJ(1) ^ 2;
t32 = sqrt(-((pkin(1) + t46) * (pkin(1) + t45) * (pkin(1) - t45) * (pkin(1) - t46)));
t1 = (((t35 * g(1) + t53 * g(2)) * t32 + ((t47 - 2 * pkin(3) * t51 + (t36 * t52 - 2 * m(2) * t42 + (-pkin(2) ^ 2 + pkin(3) ^ 2 - t42) * m(3)) * g(1)) * qJ(1)) + ((((-m(2) - m(3)) * t52 + t36) * g(1) - t51) * t42) + (pkin(3) * t47) - ((t36 * g(1) - t51) * t49) + ((-t51 + (m(3) * t52 + t53) * g(1)) * t41)) / t48 / 0.2e1 + ((rSges(4,1) * g(2) - rSges(4,2) * g(1)) * t32 - (2 * g(3) * pkin(1) * pkin(2) * rSges(4,3)) + ((qJ(1) * t52 + t41 - t42 - t49) * (rSges(4,1) * g(1) + rSges(4,2) * g(2)))) * m(4) / pkin(2) / 0.2e1) / pkin(1);
U = t1;

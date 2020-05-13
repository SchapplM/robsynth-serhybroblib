% Calculate potential energy for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
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
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fourbarprisOL_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energypot_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_energypot_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisOL_energypot_fixb_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:02
% EndTime: 2020-05-07 09:52:02
% DurationCPUTime: 0.04s
% Computational Cost: add. (21->20), mult. (30->28), div. (0->0), fcn. (14->4), ass. (0->6)
t11 = -rSges(3,3) - qJ(2);
t10 = cos(qJ(1));
t9 = cos(qJ(3));
t8 = sin(qJ(1));
t7 = sin(qJ(3));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (-t10 * rSges(2,1) + t8 * rSges(2,2) + pkin(1)) + g(2) * (-t8 * rSges(2,1) - t10 * rSges(2,2)) + g(3) * rSges(2,3)) - m(3) * (g(1) * (-t8 * rSges(3,1) + t11 * t10 + pkin(1)) + g(2) * (t10 * rSges(3,1) + t11 * t8) - g(3) * rSges(3,2)) - m(4) * (g(1) * (-t9 * rSges(4,1) + t7 * rSges(4,2)) + g(2) * (-t7 * rSges(4,1) - t9 * rSges(4,2)) + g(3) * rSges(4,3));
U = t1;

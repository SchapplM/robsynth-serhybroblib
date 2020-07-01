% Calculate minimal parameter regressor of potential energy for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% U_reg [1x9]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:21
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = fourbar1TE_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_energypot_fixb_regmin_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1TE_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:20:58
% EndTime: 2020-06-26 17:20:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (98->23), mult. (152->29), div. (8->3), fcn. (44->4), ass. (0->21)
t69 = cos(qJ(1));
t82 = pkin(2) * t69;
t81 = (-0.2e1 * t82 + pkin(1)) * pkin(1);
t79 = pkin(2) ^ 2 + t81;
t86 = 0.1e1 / t79 / 0.2e1;
t85 = -pkin(3) - pkin(4);
t84 = -pkin(3) + pkin(4);
t68 = sin(qJ(1));
t83 = pkin(2) * t68;
t80 = pkin(3) ^ 2 - pkin(4) ^ 2;
t78 = 0.1e1 / pkin(4) * t86;
t77 = 0.1e1 / pkin(3) * t86;
t76 = -pkin(1) + t82;
t65 = t79 - t80;
t64 = t79 + t80;
t63 = t76 * g(1) + g(2) * t83;
t62 = g(1) * t83 - t76 * g(2);
t61 = sqrt(-((pkin(2) - t85) * (pkin(2) + t85) + t81) * ((pkin(2) - t84) * (pkin(2) + t84) + t81));
t60 = t63 * t61;
t59 = t62 * t61;
t1 = [0, -g(1) * t69 - g(2) * t68, g(1) * t68 - g(2) * t69, 0, (t64 * t63 - t59) * t77, (-t64 * t62 - t60) * t77, 0, (-t65 * t63 - t59) * t78, (t65 * t62 - t60) * t78;];
U_reg = t1;

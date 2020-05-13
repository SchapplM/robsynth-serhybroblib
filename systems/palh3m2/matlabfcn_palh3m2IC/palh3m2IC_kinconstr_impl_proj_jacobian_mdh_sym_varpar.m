% Jacobian of implicit kinematic constraints of
% palh3m2IC
% projection from active to passive joints coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% B21 [(no of passive joints)x(no of active joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 05:00
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function B21 = palh3m2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2IC_kinconstr_impl_proj_jacobian_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_projection_jacobian_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:55:20
% EndTime: 2020-05-07 04:55:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (221->49), mult. (67->42), div. (24->8), fcn. (45->22), ass. (0->31)
t103 = (qJ(2) - qJ(6));
t100 = qJ(4) + pkin(14);
t95 = (qJ(3) + t100);
t93 = (pkin(15) + t95);
t91 = (t93 - t103);
t87 = -2 * qJ(7) - pkin(16) + t91;
t92 = (t93 + t103);
t88 = pkin(16) + t92;
t89 = -qJ(7) + t92;
t90 = -qJ(7) + t91;
t109 = (cos(t89) - cos(t90)) * pkin(1) - (cos(t87) - cos(t88)) * pkin(2);
t107 = -cos((qJ(8) - t87)) + cos((qJ(8) - t88));
t86 = 0.1e1 / pkin(2);
t108 = 0.1e1 / t107 * t86;
t102 = qJ(7) + qJ(8);
t97 = pkin(15) - t102;
t104 = cos((2 * t95)) - cos((2 * t97));
t101 = qJ(7) + pkin(16);
t94 = t101 + t103;
t78 = sin(t94);
t77 = 0.1e1 / t78;
t99 = pkin(1) * t77 / pkin(5);
t84 = 0.1e1 / pkin(7);
t98 = -pkin(4) / sin((-t93 + t102)) * t84;
t83 = 0.1e1 / pkin(9);
t96 = pkin(3) * t83 * t108;
t81 = sin(t103);
t80 = sin(t101);
t79 = sin(t100);
t66 = pkin(1) * (cos((qJ(8) - t103)) - cos((qJ(8) + t103))) + (cos((qJ(8) - t94)) - cos((qJ(8) + t94))) * pkin(2);
t1 = [0, t66 * t96, -(t104 * pkin(9) + (cos((2 * qJ(3)) + t100) - cos((2 * qJ(7)) - (2 * pkin(15)) + (2 * qJ(8)) + t100)) * pkin(4)) * t83 / t104, 0; 0, t80 * t99, 0, 0; 0, (-pkin(1) * t81 - pkin(2) * t78) * t86 * t77, 0, 0; 0, (-t109 * pkin(3) + (t107 * pkin(2) + (-cos((-qJ(8) + t90)) + cos((-qJ(8) + t89))) * pkin(1)) * pkin(7)) * t84 * t108, t79 * t98, 0; 0, (pkin(2) * t80 + pkin(5) * t81) * t86 * t99, 0, 0; 0, -(pkin(7) * t66 + t109 * pkin(9)) * t84 * t96, (pkin(9) * t79 + pkin(7) * sin((qJ(3) + t97))) * t83 * t98, 0;];
B21 = t1;
